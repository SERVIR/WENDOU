///////////////////////////////////////////////////////////////////////////////////////
///
/// Tim Mayer SCO Dev 7/19/24, prior Dev Biplov Bhandari, and Kel Markert SCO
///
/// Equations Based on Soti et al. (2030) --> https://hess.copernicus.org/articles/14/1449/2010/hess-14-1449-2010.pdf
/// Threshold Based on Herndon et al. (2018) --> https://doi.org/10.3390/s20020431
/// Workflow fo rthe methods in the script: https://docs.google.com/drawings/d/1sZj6ssqMPxiaA9j7P9dr0dtOLY2Ton8K7D_QHOuj8jY/edit
///
/// This script uses historical CHIRPS for 2019 to derive the Area-Height Relationship.
///////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////
///
/// Add FC with Ponds | Merge additional users ponds | pick pond ID of intrest
///
///////////////////////////////////////////////////////////////////////////////////////
var pondsx = ee.FeatureCollection("projects/ee-lissongdiop/assets/pithirki")

pondsx = pondsx.map(function(feature) { 
  ///////For newly user added geomtires use this function to add to the exiting ponds fc to ensure the IDs dont throw an issue
  var newDict = {Classe: 'XXX',
              Code: 999,
              Superficie: 99999,
              uniqID:360
  };
  var featureUpdate = feature.set(newDict);
  
  return featureUpdate; 
  
})

Map.addLayer(pondsx, {color:"red"}, "pondsx")
print("pondsx", pondsx)

/////////existing FC of ponds 
var pondsy = ee.FeatureCollection("users/kelmarkert/public/ferloPonds")
Map.addLayer(pondsy, {color:"green"}, "pondsy")
Map.centerObject(pondsy)
print("pondsy", pondsy)

var ponds = pondsy.merge(pondsx)
Map.addLayer(ponds, {color:"blue"}, "ponds")
print("ponds", ponds)


////////////////////
///Select the pond of intrest
var pondId = 360; /// this number needs to be updated 

var pond = ee.Feature(ponds.filter(ee.Filter.eq('uniqID', pondId)).first());
print("Sample pond geometry", pond.geometry());
Map.addLayer(pond, {color: 'red'}, 'pond');
Map.centerObject(pond, 14);

///////////////////////////////////////////////////////////////////////////////////////
///
/// Add EO
///
///////////////////////////////////////////////////////////////////////////////////////

var chirps = ee.ImageCollection("UCSB-CHG/CHIRPS/DAILY")
chirps = chirps.select(['precipitation'], ['precip']);

//var lc8 = ee.ImageCollection("LANDSAT/LC08/C02/T1_TOA") ///Removed in 7/19/24 version | if used in the future MNDWI band below need to be updated 

function maskS2clouds(image) {
  var qa = image.select('QA60');

  // Bits 10 and 11 are clouds and cirrus, respectively.
  var cloudBitMask = 1 << 10;
  var cirrusBitMask = 1 << 11;

  // Both flags should be set to zero, indicating clear conditions.
  var mask = qa.bitwiseAnd(cloudBitMask).eq(0)
      .and(qa.bitwiseAnd(cirrusBitMask).eq(0));

  return image.updateMask(mask).divide(10000).copyProperties(image, ['system:time_start']);
}

var S2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED').filterBounds(pond.geometry())
                  // .filterDate('2020-01-01', '2020-01-30')
                  // Pre-filter to get less cloudy granules.
                  // .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE',20))
                  .map(maskS2clouds);

var visualization = {
  min: 0.0,
  max: 0.3,
  bands: ['B4', 'B3', 'B2'],
};
print("S2", S2.first())

Map.addLayer(S2.mean(), visualization, 'RGB');

var elv_org = ee.ImageCollection("projects/servir-wa/SETSM_dem/SETSM_dem2").mosaic().select(["b1"], ["elevation"]);
var elv = elv_org.reproject(ee.Projection('EPSG:4326').atScale(2));
var demScale = elv.projection().nominalScale();
print("DEM resolution", demScale);

Map.addLayer(elv_org, {min: 40, max: 70}, 'dem_2m', false);


///////////////////////////////////////////////////////////////////////////////////////
///
/// Model set up | timeframe of intest and MNDWI threshold
///
///////////////////////////////////////////////////////////////////////////////////////

var forecastDays = 365;

////////update the date of interest for use
var initDate = ee.Date('2023-01-01');

var rainStudyPeriod = chirps.filterDate(initDate.advance(-20, 'day'),initDate.advance(1, 'day'));
// print('rainStudyPeriod', rainStudyPeriod);
Map.addLayer(ee.Image(rainStudyPeriod.mean()).randomVisualizer(), {}, 'rainStudyPeriod', false);

var initImg = ee.Image(S2.filterBounds(pond.geometry()).filterDate(initDate.advance(-20, 'day'),initDate.advance(1, 'day')).sort('system:time_start',false).first());
print('initImg', initImg);

var t = ee.Date(initImg.get('system:time_start'));
print('initImg time', t);

///// Awarness when swapping to Sentinel 2 MNDWI is B3 (green) and B11 (SWIR2)
///// Threshold set a -0.35
var initWater = initImg.normalizedDifference(['B3', 'B11']).gt(-0.35); // update as needed

var initPct = ee.Number(initWater.reduceRegion({
  geometry: pond.geometry(),
  reducer: ee.Reducer.mean(),
  scale: 30,
  maxPixels: 1e9
}).get('nd'));

print("initPct", initPct);

///////////////////////////////////////////////////////////////////////////////////////
///
/// Model parameters | Set Analsysis period Default to 30
///
///////////////////////////////////////////////////////////////////////////////////////
// water balance parameters [from Soti el al. (2010), Table 2]
var k = ee.Image(0.9);            // dimensionless | coefficient expressing soil moisture decrease in time | ranges: 0-1
var Gmax = ee.Image(0.01487);     // m/day | rainfall threshold value start runoff in dry soils | ranges: 0.01-0.02
var L = ee.Image(0.00114);// ***  // m/day | water loss per day | range: 0.005-0.02
var Kr = ee.Image(0.4946);        // dimensionless | runoff coefficient | range = 0.15-0.40
var n = ee.Image(19.89);          // dimensionless | # times catchment area of small pond is larger than the max pond surface area | range: 1-20
var Ac = n.multiply(pond.area()); // sq m | catchement area | range: 0-150,000,000
var alpha = ee.Image(2.514);      // dimensionless | water body shape factor | 1-3


// volume-area-height state variables
var ho = ee.Image(1); // m | pond water height | h0 state --> 1m water height

var _pondMin = elv.reduceRegion({
  geometry: pond.geometry(),
  reducer: ee.Reducer.min(),
  scale:demScale
}).get('elevation');

var pondMin = ee.Image.constant(_pondMin);
print("pond Min height is ", _pondMin.getInfo() + " meters");

var SoInit = ee.Image(0).where(elv.gte(pondMin).and(elv.lte(pondMin.add(ho))), 1);
//Map.addLayer(SoInit, {min:0, max:1}, 'SoInit', false);

var nPixels = SoInit.reduceRegion({
  geometry: pond.geometry(),
  reducer: ee.Reducer.sum(),
  scale: demScale
}).get('constant');
print('nPixels', nPixels);

var So = ee.Image(ee.Number(nPixels)).multiply(ee.Image.pixelArea());
//Map.addLayer(So, {}, 'So', false);

var SoArea = So.reduceRegion({
  geometry: pond.geometry(),
  reducer: ee.Reducer.first(),
  scale: demScale
}).get('constant');

print('The SoArea is ', SoArea.getInfo() + ' m2');

// calculate initial conditions
var A = ee.Image(pond.area()).multiply(ee.Image(initPct));               // A = pond area at time t [written as A(t)] --> this equation gives the surface area of pond actually covered in water from RS data at given time
var hInit = ho.multiply((A.divide(So)).pow(ee.Image(1).divide(alpha)));  // Eq 6 solve for h(t)
var Vo = (So.multiply(ho)).divide(alpha.add(1));                         // Eq 7 find Vo, Vo is the volume for ho=1m of water height in pond
var vInit = ee.Image(Vo.multiply((hInit.divide(ho)).pow(alpha.add(1)))); // Eq 7 find vInit based on value of Vo just calculated

var currentExtent = A.reduceRegion({geometry:pond.geometry(), reducer:ee.Reducer.first(), scale:30});
print("Maxmum extent of pond:", pond.area()," m2 and current surface area extent: ", currentExtent.get('constant'), " m2");

// set contants
// this was original scale; however this same scale works for CHIRPS coverting value from mm/day to m/day.
var precipScale = ee.Image(1).divide(ee.Image(1e3));

////

var dailyPrecip = accumChirps(chirps, t, forecastDays, precipScale);
print("Daily Precip", dailyPrecip);
Map.addLayer(dailyPrecip, {}, 'dailyPrecip', false);

///
var pastDays = 30;  ////////User update as needed

///// see the sets of Functions below as well asworkflow for better understanding of methods: https://docs.google.com/drawings/d/1sZj6ssqMPxiaA9j7P9dr0dtOLY2Ton8K7D_QHOuj8jY/edit

// InitIap is a weighted summation of past daily precip amounts used to indicate amount of water in soil
var initIap = calcInitIapWithChirps(chirps.filterDate(t.advance(-pastDays, 'day'), t), pastDays);
// print("Initial Iap", initIap);
Map.addLayer(initIap, {min:0, max:0.1}, "initIap", false);

// set initial conditions with t-1 forcing | this puts all variables into different bands of 1 image for use in volume model
var first = ee.Image(chirps.filterDate(t.advance(-1, 'day'), t).first())
              .multiply(precipScale).addBands(initIap) // converts to m/day units (but is this correct due to s^-1 of original units?)
              .addBands(vInit).addBands(A).addBands(hInit)
              .rename(['precip','Iap','vol','area','height'])//.clip(studyArea)
              .set('system:time_start', t.advance(-1, 'day').millis(), 'system:time_end', t.advance(-1, 'day').millis()).float();
print("initial variables with (t-1) forcing", first);

///////////////////////////////////////////////////////////////////////////////////////////////
//
// Initialize volume model to produce an image collection of 'precip','Iap','vol','area','height' over user specified time period through nested iterate function
//
///////////////////////////////////////////////////////////////////////////////////////////////
var modelOut = ee.ImageCollection.fromImages(dailyPrecip.iterate(accumVolume, ee.List([first])));
print("Model out", modelOut); 

Map.addLayer(modelOut.select('vol'), {}, 'modelOut', false);

// commented off (not nesscary) returns pond surface area as percentage of total area from pond shapefiles
// var pondPct = modelOut.select('area').map(function(img) {
//   var pct = img.divide(ee.Image(pond.area())).copyProperties(img, ['system:time_start']);
//   return pct;//.where(pct.gt(1),1) --> sets maximum fill at 100%
// });
// print("PondPct", pondPct);

/////
//Chart mean volume from modelOut (aka Wendou) overtime
/////
var volTimeSeries = ui.Chart.image.seriesByRegion({
  imageCollection: modelOut.limit(30), // remove limit to chart all volume model data 
  regions: pond.geometry(),
  reducer: ee.Reducer.mean(),
  band: 'vol',
  scale: demScale,
  xProperty: 'system:time_start',
  seriesProperty: 'label'
});
volTimeSeries.setChartType('ScatterChart');
volTimeSeries.setOptions({
  title: 'Pond volume vs time',
  vAxis: {
    title: 'Volume (m^3)',
  },
  lineWidth: 2,
  pointSize: 3,
  series: {
    0: 'red'
  }
});

print('Temporal Trend of the Volume', volTimeSeries);

///////////////////////////////////////////////////////////////////////////////////////////////
//
// Export model outputs for forecastDays | Produce pond specific image collection of model outputs
//
///////////////////////////////////////////////////////////////////////////////////////////////

var modelOutLists = modelOut.toList(modelOut.size());
for (var i=301; i<366; i++) {
  var img = ee.Image(modelOutLists.get(i));
  Export.image.toAsset({
    image: img,
    description: 'img_'+i,
    assetId: 'projects/servir-wa/services/ephemeral_water_ferlo/hydro_model_output/2023/pond_id_' + pondId + '_doy_' + i,
    region: pond.geometry().bounds(),
    scale: demScale,
    maxPixels: 1E13
  });
}

///////////////////////////////////////////////////////////////////////////////////////////////
//
// Wendou Functions | accumVolume | accumChirps | timeScale | accumCFS | calcInitIap | calcInitIapWithChirps | function commented out accumGFS (no longer sued past dev)
// See equations Based on Soti et al. (2030): https://hess.copernicus.org/articles/14/1449/2010/hess-14-1449-2010.pdf
// See: https://docs.google.com/drawings/d/1sZj6ssqMPxiaA9j7P9dr0dtOLY2Ton8K7D_QHOuj8jY/edit
///////////////////////////////////////////////////////////////////////////////////////////////

function accumVolume(img,list) {
  // extract out forcing and state variables
  // "past" equivalent to the x(t-1) state of variables
  var past = ee.Image(ee.List(list).get(-1));
  var pastIt = past.select('Iap');
  var pastPr = past.select('precip');
  var pastAr = past.select('area');
  var pastHt = past.select('height');
  var pastVl = past.select('vol');
  var nowPr = img.select('precip');
  var date = ee.Date(img.get('system:time_start'));

  // change in volume model
  var deltaIt = pastIt.add(pastPr).multiply(k);                                // Eq 5
  var Gt = Gmax.subtract(deltaIt);                                             // Eq 4
  Gt = Gt.where(Gt.lt(0),0);                                                   // Eq 4 (cont)
  var Pe = nowPr.subtract(Gt);                                                 // Eq 3
  Pe = Pe.where(Pe.lt(0),0);                                                   // Eq 3 (cont)
  var Qin = Kr.multiply(Pe).multiply(Ac);                                      // Eq 2
  var dV = nowPr.multiply(pond.area()).add(Qin).subtract(L.multiply(pastAr));  // Eq 1 (Qout is assumed to be 0 in Ferlo use case)

  // convert dV to actual volume (add change in volume to the initial volume to get volume at given t step)
  var volume = pastVl.add(dV).rename('vol');
  volume = volume.where(volume.lt(0), 0);

  // empirical model for volume to area/height relationship
  var ht = ho.multiply(volume.divide(Vo).pow(ee.Image(1).divide(alpha.add(1)))).rename('height');
  ht = ht.where(ht.lt(0),0);
  var area = So.multiply(ht.divide(ho).pow(alpha)).rename('area'); //Eq 6
  area = area.where(area.lt(0),1); // constrain area to real values

  // set state variables to output model step
  var step = nowPr.addBands(deltaIt).addBands(volume).addBands(area).addBands(ht)
              .set('system:time_start', date.advance(1, 'day').millis());

  return ee.List(list).add(step.float());
}

function accumChirps (collection, startDate, nDays, scale) {
  if (scale === undefined) { scale = 0 }
  // print("scale accumchirps", scale);
  // chirps has daily values in it so
  var dailyPrecip = ee.ImageCollection(collection.filterDate(startDate, startDate.advance(nDays, 'day')));
  dailyPrecip = dailyPrecip.map(function (img) {
    var sd = img.get('system:time_start');
    var ed = img.get('system:time_end');
    img = ee.Image(ee.Algorithms.If(
      scale !==0,
      img.multiply(scale).copyProperties(img).set('system:time_start', sd, 'system:time_end', ed),
      img.copyProperties(img).set('system:time_start', sd, 'system:time_end', ed)
    ));
    // return img.multiply(precipScale).copyProperties(img).set('system:time_start', sd, 'system:time_end', ed);
    return img;
  });
  return dailyPrecip;
}

function timeScale (img){
  return img.multiply(60*60*6);
}

function accumCFS(collection,s,nDays) {
  var imgList = [];
  for (var i=0; i<nDays; i++) {
    var newDate = s.advance(i,'day');
    var dayPrecip = collection.filterDate(newDate,newDate.advance(24,'hour'))
      .map(timeScale);
    imgList.push(dayPrecip.sum().multiply(precipScale)
      .set('system:time_start',s.advance(i,'day')));
  }
  return ee.ImageCollection(imgList);
}

function calcInitIap(collection, pastDays) {
  var off = pastDays*-1;
  var s = t.advance(off,'day');
  var e = s.advance(pastDays,'day');
  var prevPrecip = collection.filterDate(s,e); // these lines give you precip of select past days (ie past 7 days)

  var dailyPrev = accumCFS(prevPrecip,s,pastDays);

  var imgList = dailyPrev.toList(pastDays);
  var outList = [];  //209-220 Eq 5

  for (var i=0; i<pastDays; i++) {
    var pr = ee.Image(imgList.get(i));
    var antecedent = pr.multiply(ee.Image(1).divide(pastDays-i));
    outList.push(antecedent);
  }
  var Iap = ee.ImageCollection(outList).sum().rename('Iap');
  return Iap;
}

function calcInitIapWithChirps(collection, pastDays) {
  var off = pastDays*-1;
  var s = t.advance(off, 'day');
  var e = s.advance(pastDays, 'day');
  var prevPrecip = collection.filterDate(s, e); // these lines give you precip of select past days

  var dailyPrev = accumChirps(prevPrecip, s, pastDays, precipScale);

  var imgList = dailyPrev.toList(pastDays);
  var outList = [];

  for (var i=0; i<pastDays; i++) {
    var pr = ee.Image(imgList.get(i));
    var antecedent = pr.multiply(ee.Image(1).divide(pastDays-i));
    outList.push(antecedent);
  }
  var Iap = ee.ImageCollection(outList).sum().rename('Iap');
  return Iap;
}

// Convert to daily precip
// function accumGFS(collection,startDate,nDays) {
//   if (nDays>16){
//     alert('Max forecast days is 16, only producing forecast for 16 days...');
//     nDays = 16;
//   }
//   var cnt = 1;
//   var imgList = [];
//   for (var i=0; i<=nDays; i++) {
//     var cntMax =(24*(i+1));
//     var forecastMeta = [];
//     for(cnt;cnt<=cntMax;cnt++){forecastMeta.push(cnt)}
//     var dayPrecip = collection.filter(ee.Filter.inList('forecast_hours', forecastMeta));
//     imgList.push(dayPrecip.sum().multiply(precipScale)
//       .set('system:time_start',startDate.advance(i,'day')));
//   }
//   return ee.ImageCollection(imgList);
// }
