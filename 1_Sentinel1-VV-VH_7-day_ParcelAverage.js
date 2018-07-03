/**** Start of imports. If edited, may not auto-convert in the playground. ****/
var brp2017 = ee.FeatureCollection("users/gglemoine62/BRP_gewaspercelen_2017_concept");
/***** End of imports. If edited, may not auto-convert in the playground. *****/

///////////////////////////////////////////
// A / FUNCTIONS
///////////////////////////////////////////
// Functions to convert from/to dB
function toNatural(img) {
  return ee.Image(10.0).pow(img.select('..').divide(10.0)).copyProperties(img, ['system:time_start'])
}

function toDB(img) {
  return ee.Image(img).log10().multiply(10.0);
}

// Remove ugly edges
function maskEdge(img) {
  var mask = img.select(0).unitScale(-25, 5).multiply(255).toByte().connectedComponents(ee.Kernel.rectangle(1,1), 100);
  return img.updateMask(mask.select(0));  
}


///////////////////////////////////////////
// B / LOAD INPUTS
///////////////////////////////////////////

// 1. Date
var start_date = '2017-04-01'
var end_date = '2017-08-01'

// 2. Get Netherlands municipalities
var gemeenten = ee.FeatureCollection('ft:1B3v8wxCkO1aGd8jF4byitKEjolHvQFyMF9nZFsA8');

// 3. Set the AOI to the collections of municiplaties of interest
// (1) Gelderse Vallei area
var aoi = gemeenten.filter(ee.Filter.inList('gemnaam', ['Ede', 'Wageningen', 'Renkum', 'Barneveld', 'Arnhem', 'Putten', 'Nijkerk']));
// (2) Utrecht/Groene Hart area
// var aoi = gemeenten.filter(ee.Filter.inList('gemnaam', ['De Ronde Venen', 'Woerden', 'Breukelen', 'Maarssen', 'Soest', 'Zeist', 'Baarn', 'De Bilt']));

var step = 7 // in days (time window for meaning)


///////////////////////////////////////////
// C / CODE
///////////////////////////////////////////

brp2017 = brp2017.map(function(f) { return f.set({'id': f.id(), 'area': f.area(), 'perimeter': f.perimeter()}) })
brp2017 = brp2017.filterMetadata('area', 'less_than', 10000000) // Remove degenerate right winding polygons

// Internally buffer parcels to avoid boundary pixels
brp2017 = brp2017.map(function(f) { return f.buffer(-10)})
brp2017 = brp2017.map(function(f) { return f.set({'bufferedarea': f.area()}) })

// Clip to AOI
brp2017 = brp2017.filterBounds(aoi)

// get the data from S1 (VV pol.)
var s1 = ee.ImageCollection('COPERNICUS/S1_GRD').filterMetadata('instrumentMode', 'equals', 'IW').
  filter(ee.Filter.eq('transmitterReceiverPolarisation', ['VV', 'VH'])).
  filterBounds(aoi).filterDate(start_date, end_date).
  sort('system:time');

// Remove ugly edges
s1 = s1.map(maskEdge)

// Extracts are made from natural (non-logarithmic) values
s1 = s1.map(toNatural)

// Olha Danylo's procedure to create weekly means (adapted)

var days = ee.List.sequence(0, ee.Date(end_date).difference(ee.Date(start_date), 'day'), step).
  map(function(d) { return ee.Date(start_date).advance(d, "day") })

var dates = days.slice(0,-1).zip(days.slice(1))

var s1res = dates.map(function(range) {
  var dstamp = ee.Date(ee.List(range).get(0)).format('YYYYMMdd')
  var temp_collection = s1.filterDate(ee.List(range).get(0),
    ee.List(range).get(1)).mean().select(['VV', 'VH'], [ee.String('VV_').cat(dstamp), ee.String('VH_').cat(dstamp)])
  return temp_collection
})

// Convert ImageCollection to image stack
function stack(i1, i2)
{
  return ee.Image(i1).addBands(ee.Image(i2))
}

var s1stack = s1res.slice(1).iterate(stack, s1res.get(0))
s1stack = ee.Image(s1stack).clip(aoi)

// Export the parcel means for this image stack for use in tensorflow runs
Export.table.toDrive(ee.Image(s1stack).reduceRegions({collection: brp2017, reducer: ee.Reducer.mean(), scale: 10}).
  select(ee.List(['id', 'area', 'bufferedarea', 'perimeter', 'gws_gewasc', 'gws_gewas']).cat(ee.Image(s1stack).bandNames()), null, false), "Landsense_GH_Dump_Region1")

// Optionally, display some combination to see if all went well
Map.addLayer(toDB(ee.Image(s1res.get(0)).addBands(ee.Image(s1res.get(s1res.size().divide(2).floor()))).addBands(ee.Image(s1res.get(-1)))).clip(aoi), 
  {bands: ['VV_20170401', 'VV_20170527', 'VV_20170722'], min: -25, max: 0 }, "S1 VV stack")

Map.addLayer(toDB(ee.Image(s1res.get(0)).addBands(ee.Image(s1res.get(s1res.size().divide(2).floor()))).addBands(ee.Image(s1res.get(-1)))).clip(aoi), 
  {bands: ['VH_20170401', 'VH_20170527', 'VH_20170722'], min: -30, max: -5 }, "S1 VH stack")

