/**** Start of imports. If edited, may not auto-convert in the playground. ****/
var brp2017 = ee.FeatureCollection("users/gglemoine62/BRP_gewaspercelen_2017_concept"),
    natura2000 = ee.FeatureCollection("users/rdandrimont/Natura2000_end2016_NL");
/***** End of imports. If edited, may not auto-convert in the playground. *****/

///////////////////////////////////////////
// A / FUNCTIONS
///////////////////////////////////////////

// ADD A NDVI BAND
function addNdvi(img) {
  var nd = img.normalizedDifference(['nir', 'red']);
  return img.addBands(nd.float().rename('NDVI'));
}


function addBSI(img) {
  var bsi = img.expression('((swir1+red)-(nir+blue))/(((swir1+red)+(nir+blue))*1.0)', {
      'swir1': img.select('swir1'),
      'red': img.select('red'),
      'nir': img.select('nir'),
      'blue': img.select('blue')
}
);
    return img.addBands(bsi.float().rename('BSI'));
}


function addBare(img) {
  var bsi = img.expression('BSI>0', {
      'BSI': img.select('BSI'),
}
);
    return img.addBands(bsi.int().rename('Bare'));
}



// CONVERT SENTINEL TO TOA
function sentinel2toa(img) {
  var toa = img.select(['B1','B2','B3','B4','B6','B8','B8A','B9','B10', 'B11','B12'],  
                       ['aerosol', 'blue', 'green', 'red', 'red2','nir','red4','h2o', 'cirrus','swir1', 'swir2'])
                       .divide(10000)
                       .addBands(img.select(['QA60']))
                       .set('solar_azimuth',img.get('MEAN_SOLAR_AZIMUTH_ANGLE'))
                       .set('solar_zenith',img.get('MEAN_SOLAR_ZENITH_ANGLE'))
                       .set('date',ee.Date(img.get('system:time_start')))
                       .set('system:time_start',img.get('system:time_start'));
    return toa
}

// FLAG CLOUD AND CIRRUS
function ESAcloud(toa) {
  // author: Nick Clinton
  var qa = toa.select('QA60');
  
  // Bits 10 and 11 are clouds and cirrus, respectively.
  var cloudBitMask = Math.pow(2, 10);
  var cirrusBitMask = Math.pow(2, 11);
  
  // clear if both flags set to zero.
  var clear = qa.bitwiseAnd(cloudBitMask).eq(0).and(
             qa.bitwiseAnd(cirrusBitMask).eq(0));
    var cloud = clear.eq(0)
  return cloud
}

// FLAG SHADOW
function shadowMask(toa,cloud){
  // Author: Gennadii Donchyts
  // License: Apache 2.0
  
  // solar geometry (radians)
  var azimuth =ee.Number(toa.get('solar_azimuth')).multiply(Math.PI).divide(180.0).add(ee.Number(0.5).multiply(Math.PI));
  var zenith  =ee.Number(0.5).multiply(Math.PI ).subtract(ee.Number(toa.get('solar_zenith')).multiply(Math.PI).divide(180.0));

  // find where cloud shadows should be based on solar geometry
  var nominalScale = cloud.projection().nominalScale();
  var cloudHeights = ee.List.sequence(200,10000,500);
  var shadows = cloudHeights.map(function(cloudHeight){
    cloudHeight = ee.Number(cloudHeight);
    var shadowVector = zenith.tan().multiply(cloudHeight);
    var x = azimuth.cos().multiply(shadowVector).divide(nominalScale).round();
    var y = azimuth.sin().multiply(shadowVector).divide(nominalScale).round();
    return cloud.changeProj(cloud.projection(), cloud.projection().translate(x, y));
  });
  var potentialShadow = ee.ImageCollection.fromImages(shadows).max();
  
  // shadows are not clouds
  var potentialShadow = potentialShadow.and(cloud.not());
  
  // (modified by Sam Murphy) dark pixel detection 
  var darkPixels = toa.normalizedDifference(['green', 'swir2']).gt(0.25).rename(['dark_pixels']);
  
  // shadows are dark
  var shadow = potentialShadow.and(darkPixels).rename('shadows');
  
  return shadow
}

// CONVERT TO TOA, MASK SHADOW AND MASK CLOUD
function cloud_and_shadow_mask(img) {
    var toa = sentinel2toa(img)
    var cloud = ESAcloud(toa)
    var shadow = shadowMask(toa,cloud)
    var mask = cloud.or(shadow).eq(0)
    return toa.updateMask(mask)
  }


// This function adds a band representing the image timestamp.
var addTime = function(image) {
  return image.addBands(image.metadata('system:time_start')
    .divide(1000 * 60 * 60 * 24 * 365));
};

///////////////////////////////////////////
// B / LOAD INPUTS
///////////////////////////////////////////

// 1. Date
var startDateStr='2017-01-01'
var stopDateStr='2017-08-01'

// 2. Get Netherlands municipalities
var gemeenten = ee.FeatureCollection('ft:1B3v8wxCkO1aGd8jF4byitKEjolHvQFyMF9nZFsA8');

// 3. Set the AOI to the collections of municiplaties of interest
// (1) Gelderse Vallei area
var aoi = gemeenten.filter(ee.Filter.inList('gemnaam', ['Ede', 'Wageningen', 'Renkum', 'Barneveld', 'Arnhem', 'Putten', 'Nijkerk']));
// (2) Utrecht/Groene Hart area
// var aoi = gemeenten.filter(ee.Filter.inList('gemnaam', ['De Ronde Venen', 'Woerden', 'Breukelen', 'Maarssen', 'Soest', 'Zeist', 'Baarn', 'De Bilt']));


///////////////////////////////////////////
// C / CODE
///////////////////////////////////////////

// C.0 LOAD PARCELS
brp2017 = brp2017.map(function(f) { return f.set({'id': f.id(), 'area': f.area(), 'perimeter': f.perimeter()}) })
brp2017 = brp2017.filterMetadata('area', 'less_than', 10000000) // Remove degenerate right winding polygons

// Internally buffer parcels to avoid boundary pixels
brp2017 = brp2017.map(function(f) { return f.buffer(-10)})
brp2017 = brp2017.map(function(f) { return f.set({'bufferedarea': f.area()}) })

// Clip to AOI
brp2017 = brp2017.filterBounds(aoi)


// C.1  COMPOSITE
///////////////////////////////////////////
// LOAD Sentinel-2 collections for a given temporal window
var startDate = ee.Date(startDateStr)
var stopDate = ee.Date(stopDateStr)
var images = ee.ImageCollection('COPERNICUS/S2')
  .filterDate(startDate, stopDate).filterBounds(aoi);
  // .limit(10);
  
// CONVERT TO TOA AND APPLY THE SHADOW AND CLOUD MASK
var s2_cleaned = images.map(cloud_and_shadow_mask);

// C.2 ADD BANDS
///////////////////////////////////////////

// ADD NDVI BANDS
var s2_cleaned = s2_cleaned.map(addNdvi);

// ADD BSI BAND
var s2_cleaned = s2_cleaned.map(addBSI);

// ADD Bare BAND
var s2_cleaned = s2_cleaned.map(addBare);

// C.3 BSI SUM
///////////////////////////////////////////
// SUM OF BARE EVENT
var result = s2_cleaned.select('Bare').reduce('sum') ;

// C.4 / LINEAR REGRESSION
///////////////////////////////////////////

// USE NDVI and TIME
var s2_ts = s2_cleaned.map(
  function(image) {
    // Rename that band to something appropriate
    return  image.select( ['NDVI']).set('system:time_start', image.get('system:time_start'));
  }
);

var s2_ts = s2_ts.map(addTime);

// Compute the linear trend over time.
var trend = s2_ts.select(['system:time_start', 'NDVI']).reduce(ee.Reducer.linearFit());


// C.5 /Combine all results
///////////////////////////////////////////
var s2stack = result.addBands(trend.select('scale')).addBands(trend.select('offset'))


// C.6 Export the parcel means for this image stack for use in tensorflow runs
///////////////////////////////////////////
Export.table.toDrive(ee.Image(s2stack).reduceRegions({collection: brp2017, reducer: ee.Reducer.mean(), scale: 10}).
  select(ee.List(['id', 'area', 'bufferedarea', 'perimeter', 'gws_gewasc', 'gws_gewas']).cat(ee.Image(s2stack).bandNames()), null, false), "Landsense_S2_Dump_Region1")





