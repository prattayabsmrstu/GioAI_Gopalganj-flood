// =============================================================================
// SAR Flood Frequency Mapping for Gopalganj, Bangladesh (FIXED)
// Google Earth Engine JavaScript API
// Produces: Annual flood frequency maps (2019-2024) + multi-year composite
// No ground truth required | Self-supervised anomaly detection | GEE-only core
// =============================================================================

// ---------------------------------------------------------
// 1. STUDY AREA: Gopalganj District (Admin Level 2)
// ---------------------------------------------------------
var gopalganj = ee.FeatureCollection('FAO/GAUL/2015/level2')
  .filter(ee.Filter.eq('ADM1_NAME', 'Dhaka'))
  .filter(ee.Filter.eq('ADM2_NAME', 'Gopalganj'))
  .first()
  .geometry();

var roi = gopalganj.buffer(5000); // 5km buffer for edge effects
Map.centerObject(roi, 10);
Map.addLayer(roi, {color: 'yellow'}, 'ROI: Gopalganj + 5km buffer');

// ---------------------------------------------------------
// 2. SENTINEL-1 GRD COLLECTION (VH POLARIZATION, IW MODE)
// ---------------------------------------------------------
var s1 = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterBounds(roi)
  .filterDate('2019-01-01', '2025-01-01')  // 2019-2024 inclusive
  .filter(ee.Filter.eq('instrumentMode', 'IW'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
  .filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING')) // Consistent geometry
  .select('VH');

// Check image count
print('Total S1 VH images (2019-2024):', s1.size());

// ---------------------------------------------------------
// 3. SPECKLE REDUCTION: Focal Mean (3x3 kernel = 30m)
// ---------------------------------------------------------
function reduceSpeckle(img) {
  return img.focal_mean(30, 'square', 'meters')  // 30m = 3x3 kernel on 10m data
    .copyProperties(img, ['system:time_start', 'system:index', 'orbitProperties_pass']);
}

var s1_smooth = s1.map(reduceSpeckle);

// ---------------------------------------------------------
// 4. DRY-SEASON BASELINE: Median + MAD (Nov–Apr, non-monsoon)
// FIXED: calendarRange doesn't wrap year boundary, use explicit month list
// ---------------------------------------------------------
var dryMonths = [11, 12, 1, 2, 3, 4];
var dryFilter = ee.Filter.or.apply(null, dryMonths.map(function(m) {
  return ee.Filter.calendarRange(m, m, 'month');
}));

var dryCollection = s1_smooth.filter(dryFilter);

// Compute median (baseline backscatter) and MAD (robust spread)
var dryMedian = dryCollection.median().rename('VH_dry_median');
var dryMAD = dryCollection.map(function(img) {
  return img.subtract(dryMedian).abs();
}).median().multiply(1.4826).rename('VH_dry_MAD'); // 1.4826 = Gaussian consistency constant

// Visualize baseline
Map.addLayer(dryMedian.clip(roi), {min: -25, max: -5, palette: ['blue', 'white', 'red']}, 'Dry Median VH');
Map.addLayer(dryMAD.clip(roi), {min: 0, max: 5, palette: ['green', 'yellow', 'red']}, 'Dry MAD VH');

print('Dry-season images used:', dryCollection.size());

// ---------------------------------------------------------
// 5. PER-IMAGE ANOMALY DETECTION (Self-Supervised)
// ---------------------------------------------------------
// Water = significant drop in VH backscatter: (VH - median) / MAD < threshold
var ANOMALY_THRESHOLD = -2.5; // Tune this with DFO validation

function detectWater(img) {
  var anomaly = img.subtract(dryMedian).divide(dryMAD).rename('anomaly');
  var water = anomaly.lt(ANOMALY_THRESHOLD).rename('water');
  return img.addBands(anomaly).addBands(water)
    .copyProperties(img, ['system:time_start', 'system:index']);
}

var s1_anomaly = s1_smooth.map(detectWater);

// ---------------------------------------------------------
// 6. ANNUAL FLOOD FREQUENCY (per year 2019-2024)
// ---------------------------------------------------------
var years = ee.List.sequence(2019, 2024);

var annualFrequency = ee.ImageCollection(years.map(function(year) {
  var start = ee.Date.fromYMD(year, 1, 1);
  var end = start.advance(1, 'year');

  var yearly = s1_anomaly.filterDate(start, end);
  var count = yearly.select('water').sum().rename('water_count');
  var total = yearly.select('water').count().rename('total_obs');

  var frequency = count.divide(total).multiply(100).rename('flood_frequency_pct');

  return frequency.clip(roi).set({  // CLIP TO ROI
    'year': year,
    'system:time_start': start.millis()
  });
}));

print('Annual frequency collection:', annualFrequency);

// ---------------------------------------------------------
// 7. MULTI-YEAR COMPOSITE FLOOD FREQUENCY (2019-2024)
// ---------------------------------------------------------
var allWater = s1_anomaly.select('water');
var totalCount = allWater.sum().rename('total_water_count');
var totalObs = allWater.count().rename('total_observations');

var multiYearFreq = totalCount.divide(totalObs).multiply(100).rename('flood_frequency_2019_2024')
  .clip(roi);

// Visualization
var freqVis = {
  min: 0, max: 100,
  palette: ['ffffff', 'ffffb2', 'fecc5c', 'fd8d3c', 'f03b20', 'bd0026']
};

Map.addLayer(multiYearFreq, freqVis, 'Flood Frequency 2019-2024 (%)');

// ---------------------------------------------------------
// 8. EXPORT TO GOOGLE DRIVE
// ---------------------------------------------------------
// Multi-year composite
Export.image.toDrive({
  image: multiYearFreq.toFloat(),
  description: 'Gopalganj_FloodFrequency_2019_2024',
  folder: 'GEE_Flood_Gopalganj',
  fileNamePrefix: 'Gopalganj_FloodFreq_2019_2024',
  region: roi,
  scale: 20,
  crs: 'EPSG:4326',
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});

// Annual exports — FIXED: use evaluate() for client-side loop
annualFrequency.evaluate(function(col) {
  col.features.forEach(function(feat) {
    var year = feat.properties.year;
    var img = ee.Image(annualFrequency.filter(ee.Filter.eq('year', year)).first());
    Export.image.toDrive({
      image: img.toFloat(),
      description: 'Gopalganj_FloodFreq_' + year,
      folder: 'GEE_Flood_Gopalganj',
      fileNamePrefix: 'Gopalganj_FloodFreq_' + year,
      region: roi,
      scale: 20,
      crs: 'EPSG:4326',
      maxPixels: 1e13,
      fileFormat: 'GeoTIFF'
    });
  });
});

// ---------------------------------------------------------
// 9. VALIDATION PREP: Export water occurrence for DFO comparison
// ---------------------------------------------------------
Export.image.toDrive({
  image: totalCount.toInt16(),
  description: 'Gopalganj_WaterCount_2019_2024',
  folder: 'GEE_Flood_Gopalganj',
  fileNamePrefix: 'Gopalganj_WaterCount_2019_2024',
  region: roi,
  scale: 20,
  crs: 'EPSG:4326',
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});

Export.image.toDrive({
  image: totalObs.toInt16(),
  description: 'Gopalganj_TotalObs_2019_2024',
  folder: 'GEE_Flood_Gopalganj',
  fileNamePrefix: 'Gopalganj_TotalObs_2019_2024',
  region: roi,
  scale: 20,
  crs: 'EPSG:4326',
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});

// ---------------------------------------------------------
// 10. SUMMARY STATISTICS (for paper tables)
// ---------------------------------------------------------
var stats = multiYearFreq.reduceRegion({
  reducer: ee.Reducer.mean().combine(ee.Reducer.stdDev(), '', true)
    .combine(ee.Reducer.minMax(), '', true)
    .combine(ee.Reducer.percentile([25, 50, 75, 90, 95]), '', true),
  geometry: roi,
  scale: 100, // coarser for speed
  maxPixels: 1e9
});

print('Flood Frequency Stats (2019-2024):', stats);

// Area calculation (sq km)
var floodArea = multiYearFreq.gt(10).multiply(ee.Image.pixelArea()).divide(1e6)
  .reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: roi,
    scale: 20,
    maxPixels: 1e13
  });

print('Area with >10% flood frequency (sq km):', floodArea);

// ---------------------------------------------------------
// 11. DFO VALIDATION (Approximate events for Madhumati basin)
// ---------------------------------------------------------
var gopalganjFeature = ee.Feature(
  ee.FeatureCollection('FAO/GAUL/2015/level2')
    .filter(ee.Filter.eq('ADM1_NAME', 'Dhaka'))
    .filter(ee.Filter.eq('ADM2_NAME', 'Gopalganj'))
    .first()
);

Map.addLayer(gopalganjFeature, {color: 'red', fillColor: '00000000'}, 'Gopalganj boundary');

// Approximate DFO monsoon flood extents (replace with real DFO polygons later)
var dfoApprox = ee.FeatureCollection([
  ee.Feature(ee.Geometry.Rectangle([89.62, 22.82, 90.08, 23.18]), {
    event: '2019_monsoon', start: '2019-07-10', end: '2019-07-13'
  }),
  ee.Feature(ee.Geometry.Rectangle([89.58, 22.78, 90.12, 23.22]), {
    event: '2020_monsoon', start: '2020-06-30', end: '2020-07-30'
  }),
  ee.Feature(ee.Geometry.Rectangle([89.66, 22.86, 90.04, 23.14]), {
    event: '2022_monsoon', start: '2022-06-10', end: '2022-06-21'
  })
]).map(function (f) {
  return ee.Feature(f.geometry().intersection(gopalganjFeature.geometry()))
    .copyProperties(f, ['event', 'start', 'end']);
});

Map.addLayer(dfoApprox, {color: '00FFFF'}, 'Approximate DFO events');

// Per-event statistics against the flood frequency map
var eventStats = dfoApprox.map(function (f) {
  var geom = f.geometry();

  var meanFreq = multiYearFreq.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: geom,
    scale: 20,
    maxPixels: 1e9
  }).get('flood_frequency_2019_2024');

  var areaHigh = multiYearFreq.gte(50)
    .multiply(ee.Image.pixelArea()).divide(1e6)
    .reduceRegion({
      reducer: ee.Reducer.sum(),
      geometry: geom,
      scale: 20,
      maxPixels: 1e13
    }).get('flood_frequency_2019_2024');

  var areaAll = ee.Image.pixelArea().divide(1e6)
    .reduceRegion({
      reducer: ee.Reducer.sum(),
      geometry: geom,
      scale: 20,
      maxPixels: 1e13
    }).get('area');

  return f.set({
    'mean_flood_freq_pct': meanFreq,
    'area_gt50_km2': areaHigh,
    'event_area_km2': areaAll
  });
});

print('Event overlap statistics:', eventStats);

Export.table.toDrive({
  collection: eventStats,
  description: 'Gopalganj_DFO_overlap_approx',
  folder: 'GEE_Flood_Gopalganj',
  fileNamePrefix: 'Gopalganj_DFO_overlap_approx',
  fileFormat: 'CSV'
});

// ---------------------------------------------------------
// 12. AUXILIARY EXPORTS (Dry baseline, boundary, DFO events)
// ---------------------------------------------------------
// Dry Median Export
Export.image.toDrive({
  image: dryMedian.clip(roi).toFloat(),
  description: 'Gopalganj_DryMedian_VH',
  folder: 'GEE_Flood_Gopalganj',
  region: roi,
  scale: 20,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});

// Dry MAD Export
Export.image.toDrive({
  image: dryMAD.clip(roi).toFloat(),
  description: 'Gopalganj_DryMAD_VH',
  folder: 'GEE_Flood_Gopalganj',
  region: roi,
  scale: 20,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});

// Boundary Export (KML)
Export.table.toDrive({
  collection: ee.FeatureCollection([ee.Feature(roi)]),
  description: 'Gopalganj_Boundary',
  folder: 'GEE_Flood_Gopalganj',
  fileFormat: 'KML'
});

// DFO Events Export (GeoJSON)
Export.table.toDrive({
  collection: dfoApprox,
  description: 'Gopalganj_DFO_Events',
  folder: 'GEE_Flood_Gopalganj',
  fileFormat: 'GeoJSON'
});

// ---------------------------------------------------------
// END OF SCRIPT
// ---------------------------------------------------------
// Run in GEE Code Editor: https://code.earthengine.google.com/
// 1. Paste this entire script
// 2. Click Run
// 3. Check Tasks tab → Run exports
// 4. Download GeoTIFFs from Google Drive folder "GEE_Flood_Gopalganj"
// 5. Use validate.ipynb for DFO validation & publication figures