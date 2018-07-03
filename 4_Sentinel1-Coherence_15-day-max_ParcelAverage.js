/**** Start of imports. If edited, may not auto-convert in the playground. ****/
var brp2017 = ee.FeatureCollection("users/gglemoine62/BRP_gewaspercelen_2017_concept"),
    s1coh = ee.ImageCollection("users/gglemoine62/nl_coh");
/***** End of imports. If edited, may not auto-convert in the playground. *****/
brp2017 = brp2017.map(function(f) { return f.set({'id': f.id(), 'area': f.area(), 'perimeter': f.perimeter()}) })

s1coh = s1coh.map(function(f) {
  return f.set('system:time_start', ee.Date.parse('YYYYMMdd', ee.String(f.id()).split('_').get(3)))
})

brp2017 = brp2017.filterMetadata('area', 'less_than', 10000000) // Remove degenerate right winding polygons
brp2017 = brp2017.map(function(f) { return f.buffer(-10)})
brp2017 = brp2017.map(function(f) { return f.set({'bufferedarea': f.area()}) })

var gemeenten = ee.FeatureCollection('ft:1B3v8wxCkO1aGd8jF4byitKEjolHvQFyMF9nZFsA8');

var aoi = gemeenten.filter(ee.Filter.inList('gemnaam', ['Ede', 'Wageningen', 'Renkum', 'Barneveld', 'Arnhem', 'Putten', 'Nijkerk']));
// var aoi = gemeenten.filter(ee.Filter.inList('gemnaam', ['De Ronde Venen', 'Woerden', 'Breukelen', 'Maarssen', 'Soest', 'Zeist', 'Baarn', 'De Bilt']));

brp2017 = brp2017.filterBounds(aoi)

Map.centerObject(aoi, 12);

var start_date = '2017-01-01'
var end_date = '2017-08-01'


// Olha's idea to create mean images
var step = 15 // in days

var days = ee.List.sequence(0, ee.Date(end_date).difference(ee.Date(start_date), 'day'), step).
  map(function(d) { return ee.Date(start_date).advance(d, "day") })

var dates = days.slice(0,-1).zip(days.slice(1))

var s1coh_res = dates.map(function(range) {
  var dstamp = ee.Date(ee.List(range).get(0)).format('YYYYMMdd')
  var temp_collection = s1coh.filterDate(ee.List(range).get(0),
    ee.List(range).get(1)).max()
  return temp_collection.select(temp_collection.bandNames(), [ee.String('COHVV_').cat(dstamp), ee.String('COHVH_').cat(dstamp)])
  
})

function stack(i1, i2)
{
  return ee.Image(i1).addBands(ee.Image(i2))
}

var s1stack = s1coh_res.slice(1).iterate(stack, s1coh_res.get(0))

Map.addLayer(ee.Image(s1stack).select([26,12,0]).clip(aoi), {max: 0.8}, "coherence")

Export.table.toDrive(ee.Image(s1stack).reduceRegions({collection: brp2017, reducer: ee.Reducer.mean(), scale: 10}).
  select(ee.List(['id', 'area', 'bufferedarea', 'perimeter', 'gws_gewasc', 'gws_gewas']).cat(ee.Image(s1stack).bandNames()), null, false), "Landsense_GH_COH_Dump_Region1")

