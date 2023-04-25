/**** Start of imports. If edited, may not auto-convert in the playground. ****/
var ponds = ee.FeatureCollection("users/kelmarkert/public/ferloPonds"),
    mk_pond = /* color: #d63000 */ee.Feature(
        ee.Geometry.Polygon(
            [[[103.11428283212649, 16.181837537330324],
              [103.11576341150271, 16.181837537330324],
              [103.11509822366702, 16.183238846689054],
              [103.11411117074954, 16.182929735185496]]]),
        {
          "system:index": "0"
        }),
    studyArea = 
    /* color: #d63000 */
    /* shown: false */
    ee.Geometry.Polygon(
        [[[-15.866, 16.49],
          [-15.866, 14.193],
          [-12.99, 14.193],
          [-12.99, 16.49]]]),
    chirps = ee.ImageCollection("UCSB-CHG/CHIRPS/DAILY"),
    volumne_pt = /* color: #d63000 */ee.Geometry.MultiPoint(),
    lc8 = ee.ImageCollection("LANDSAT/LC08/C02/T1_TOA"),
    wendou = ee.ImageCollection("users/biplovbhandari/UAH/Wendou_2019");
/***** End of imports. If edited, may not auto-convert in the playground. *****/
// Original author: K. Markert
// Based on Soti et al. (2010) --> https://hess.copernicus.org/articles/14/1449/2010/hess-14-1449-2010.pdf
// Edited by: Biplov Bhandari . SCO (4/25/2023)
// This script uses historical CHIRPS for 2019 to derive the Area-Height Relationship.

var elv_org = ee.ImageCollection("projects/servir-wa/SETSM_dem/SETSM_dem2").mosaic().select(["b1"], ["elevation"]);
var elv = elv_org.reproject(ee.Projection('EPSG:4326').atScale(2));
var demScale = elv.projection().nominalScale(); 
print("DEM resolution", demScale);


Map.addLayer(elv, {min: 40, max: 70}, 'dem_2m');
// var studyArea = ee.Geometry.Rectangle([-180,-60,180,85])//mk_pond.buffer(10000,100).geometry()

var forecastDays = 365;
var initDate = ee.Date('2021-01-01');//ee.Date(date); 

var pondId = 75; /// test case was #1

var pond = ee.Feature(ponds.filter(ee.Filter.eq('uniqID', pondId)).first());
print("Sample pond geometry", pond.geometry());
Map.centerObject(pond, 14);

var initImg = ee.Image(lc8.filterBounds(pond.geometry()).filterDate(initDate.advance(-20, 'day'),initDate.advance(1, 'day')).sort('system:time_start',false).first());
print('initImg', initImg);
var t = ee.Date(initImg.get('system:time_start'));
print('initImg time', t);

// MNDWI on L8 with B3 (green) and B7 (SWIR2)
var initWater = initImg.normalizedDifference(['B3', 'B7']).gt(-0.2); // >-0.2 is a low threshold for MNDWI, potential overestimate of water?
Map.addLayer(initImg.normalizedDifference(['B3', 'B7']), {bands: 'nd', min: -0.4, max: -0.01}, 'MNDWI', false);

var initPct = ee.Number(initWater.reduceRegion({
  geometry: pond.geometry(),
  reducer: ee.Reducer.mean(),
  scale: 30,
  maxPixels: 1e9
}).get('nd'));

Map.addLayer(initImg.select(['B4','B3','B2']),{min:0,max:0.3,gamma:1.3},'Natural-color',false);
Map.addLayer(initImg.normalizedDifference(['B3','B7']).rename("MNDWI"),{min:-0.4,max:-0.1},'Initial Conditions',false);

print("initPct", initPct);

/* ----- model parameterization ----- */
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
Map.addLayer(SoInit, {min:0, max:1}, 'SoInit', false);

var nPixels = SoInit.reduceRegion({
  geometry: pond.geometry(),
  reducer: ee.Reducer.sum(),
  scale: demScale
}).get('constant');
print('nPixels', nPixels);

var So = ee.Image(ee.Number(nPixels)).multiply(ee.Image.pixelArea());
// print("So",So);

var SoArea = So.reduceRegion({
  geometry: pond.geometry(),
  reducer: ee.Reducer.first(),
  scale: demScale
}).get('constant');

print('The SoArea is ', SoArea.getInfo() + ' m2');

// calculate initial conditions
var A = ee.Image(pond.area()).multiply(ee.Image(initPct));               // A = pond area at time t [written as A(t)] --> this equation gives the surface area of pond actually covered in water from RS data at given time
var hInit = ho.multiply(A.divide(So).pow(ee.Image(1).divide(alpha)));       // Eq 6 solve for h(t)
var Vo = (So.multiply(ho)).divide(alpha.add(1));                         // Eq 7 find Vo, Vo is the volume for ho=1m of water height in pond
var vInit = ee.Image((Vo.multiply(hInit.divide(ho)).pow(alpha.add(1)))); // Eq 7 find vInit based on value of Vo just calculated

var currentExtent = A.reduceRegion({geometry:pond.geometry(), reducer:ee.Reducer.first(), scale:30});
print("Maxmum extent of pond:", pond.area()," m2 and current surface area extent: ", currentExtent.get('constant'), " m2");

// set contants 
// this was original scale; however this same scale works for CHIRPS coverting value from mm/day to m/day.
var precipScale = ee.Image(1).divide(ee.Image(1e3));

/* ----- start proccessing ----- */

chirps = chirps.select(['precipitation'], ['precip']);
var precipData = chirps.filterDate(t, t.advance(1, 'day')).filterBounds(studyArea);
print("Precipitation Data", precipData);
             
// var dailyPrecip = accumChirps(precipData, t, forecastDays);
var dailyPrecip = accumChirps(chirps, t, forecastDays);
dailyPrecip = dailyPrecip.map(function (img) {
  var sd = img.get('system:time_start');
  var ed = img.get('system:time_end');
  return img.multiply(precipScale).copyProperties(img).set('system:time_start', sd, 'system:time_end', ed);
});
print("Daily Precip", dailyPrecip);

var pastDays = 7;
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

///////////////////////////////////////////////////////////////////////////////////////////////  // Initialize volume model              
var modelOut = ee.ImageCollection.fromImages(dailyPrecip.iterate(accumVolume, ee.List([first])));
print("Model out", modelOut); // why is res of output the res of GFS (~27km2), rather than 30m Landsat res??

// returns pond surface area as percentage of total area from pond shapefiles
var pondPct = modelOut.select('area').map(function(img) {
  var pct = img.divide(ee.Image(pond.area())).copyProperties(img, ['system:time_start']);
  return pct;//.where(pct.gt(1),1) --> sets maximum fill at 100%
});
// print("PondPct", pondPct);


// var modelOutLists = modelOut.toList(modelOut.size());

// for (var i=300; i<=forecastDays; i++) {
//   var img = ee.Image(modelOutLists.get(i));
//   Export.image.toAsset({
//     image: img,
//     description: 'img_'+i,
//     assetId: 'users/biplovbhandari/UAH/Wendou_2019/image_' + i,
//     region: pond.geometry().bounds(),
//     scale: demScale,
//     maxPixels: 1E13
//   });
// }

print('Temporal Trend of the Surface Area Extent');
var timeSeries = ui.Chart.image.seriesByRegion({
  imageCollection: wendou,
  regions: pond.geometry(),
  reducer: ee.Reducer.mean(),
  band: 'area',
  scale: demScale,
  xProperty: 'system:time_start',
  seriesProperty: 'label'
});
timeSeries.setChartType('ScatterChart');
timeSeries.setOptions({
  title: 'Pond Area vs time',
  vAxis: {
    title: 'Area (m^2)',
  },
  lineWidth: 2,
  pointSize: 3,
  series: {
    0: 'red'
  }
});

print(timeSeries);



wendou = wendou.map(function (img) {
  var h = img.select('height');
  var log_H = h.log10().rename('log_H');
  var a = img.select('area');
  var log_A = a.log10().rename('log_A');
  var vol = img.select('vol');
  var log_V = vol.log10().rename('log_V');
  return img.addBands(log_H).addBands(log_A).addBands(log_V).addBands(ee.Image.constant(0.001).rename('log_C'));
});

// Areal Trend Line

// trend line would be
// A = C (h) ^ alpha
// this would be reduced to
// log A = alpha * log h + log C
var independents = ee.List(['log_C', 'log_H']);
var dependent = ee.String('log_A');

// Compute a linear trend.  This will have two bands: 'residuals' and 
// a 2x1 (Array Image) band called 'coefficients'.
// (Columns are for dependent variables)
var trend = wendou.select(independents.add(dependent))
    .reduce(ee.Reducer.linearRegression(independents.length(), 1));

// Flatten the coefficients into a 2-band image.
var coefficients = trend.select('coefficients')
    // Get rid of extra dimensions and convert back to a regular image
    .arrayProject([0])
    .arrayFlatten([independents]);

// Compute fitted values.
wendou = wendou.map(function(image) {
    return image.addBands(
        image.select(independents)
        .multiply(coefficients)
        .reduce('sum')
        .rename('log_fitted_area'));
});

wendou = wendou.map(function (img) {
  return img.addBands(ee.Image(10).pow(img.select('log_fitted_area')).rename('area_modeled'));
});


// Volumetric Trend Line
// trend line would be
// V = C (h) ^ alpha
// this would be reduced to
// log V = alpha * log h + log C
var independents = ee.List(['log_C', 'log_H']);
var dependent = ee.String('log_V');

// Compute a linear trend.  This will have two bands: 'residuals' and 
// a 2x1 (Array Image) band called 'coefficients'.
// (Columns are for dependent variables)
var trend = wendou.select(independents.add(dependent))
    .reduce(ee.Reducer.linearRegression(independents.length(), 1));

// Flatten the coefficients into a 2-band image.
var coefficients = trend.select('coefficients')
    // Get rid of extra dimensions and convert back to a regular image
    .arrayProject([0])
    .arrayFlatten([independents]);

// Compute fitted values.
wendou = wendou.map(function(image) {
    return image.addBands(
        image.select(independents)
        .multiply(coefficients)
        .reduce('sum')
        .rename('log_fitted_vol'));
});

wendou = wendou.map(function (img) {
  return img.addBands(ee.Image(10).pow(img.select('log_fitted_vol')).rename('vol_modeled'));
});


print('wendou', wendou);
    

var wendouList = wendou.toList(wendou.size());
var loc = pond.geometry().bounds().centroid(1);
var samples = wendouList.map(function (img) {
  img = ee.Image(img);
  return img.sample({
    region: loc,
    scale: demScale
  });
});
samples = ee.FeatureCollection(samples.flatten());
samples = samples.flatten();


// Prepare the chart.
var VAHChart =
  ui.Chart.feature.groups(samples, 'height', 'area', 'series')
    .setChartType('ScatterChart')
    .setOptions({
      title: 'Pond: ID ' + pondId,
      hAxis: {
        title: 'Height'
      },
      vAxis: {
        title: 'Area'
      },
      pointSize: 3,
    // trendlines: {
    //     0: {
    //       type: 'exponential',
    //       visibleInLegend: true,
    //       color: 'red',
    //       lineWidth: 7,
    //       opacity: 0.4,
    //     }
    // }
  });

print('Area-height curve', VAHChart);

// Prepare the chart.
var VAHChart1 =
  ui.Chart.feature.groups(samples, 'area', 'area_modeled', 'series')
    .setChartType('ScatterChart')
    .setOptions({
      title: 'Pond: ID ' + pondId,
      hAxis: {
        title: 'Area'
      },
      vAxis: {
        title: 'Fitted'
      },
      pointSize: 3,
      trendlines: {
            0: {
                color: 'red',
                lineWidth: 7,
                opacity: 0.5,
            }
        },
    });

print('Area vs fitted', VAHChart1);

var areaDiff = function(feature) {
  var diff = ee.Number(feature.get('area')).subtract(ee.Number(feature.get('area_modeled')));
  // Return the feature with the squared difference set to the 'diff' property.
  return feature.set('area_diff', diff.pow(2));
};

var rmse = ee.Number(
  // Map the difference function over the collection.
  samples.map(areaDiff)
  // Reduce to get the mean squared difference.
  .reduceColumns(ee.Reducer.mean(), ['area_diff'])
  .get('mean')
)
// Compute the square root of the mean square to get RMSE.
.sqrt();

// Print the result.
print('Area RMSE=', rmse);

// volume-height
// Prepare the chart.
var VAHChart =
  ui.Chart.feature.groups(samples, 'height', 'vol', 'series')
    .setChartType('ScatterChart')
    .setOptions({
      title: 'Pond: ID ' + pondId,
      hAxis: {
        title: 'Height'
      },
      vAxis: {
        title: 'Volume'
      },
      pointSize: 3,
    // trendlines: {
    //     0: {
    //       type: 'exponential',
    //       visibleInLegend: true,
    //       color: 'red',
    //       lineWidth: 7,
    //       opacity: 0.4,
    //     }
    // }
  });



print('Temporal Trend of the Volume');

var timeSeries = ui.Chart.image.seriesByRegion({
  imageCollection: wendou,
  regions: pond.geometry(),
  reducer: ee.Reducer.mean(),
  band: 'vol',
  scale: demScale,
  xProperty: 'system:time_start',
  seriesProperty: 'label'
});
timeSeries.setChartType('ScatterChart');
timeSeries.setOptions({
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

print(timeSeries);

print('Volume-Height Chart', VAHChart);

// Prepare the chart.
var VAHChart1 =
  ui.Chart.feature.groups(samples, 'vol', 'vol_modeled', 'series')
    .setChartType('ScatterChart')
    .setOptions({
      title: 'Pond: ID ' + pondId,
      hAxis: {
        title: 'Volume'
      },
      vAxis: {
        title: 'Fitted Volume'
      },
      pointSize: 3,
      trendlines: {
            0: {
                color: 'red',
                lineWidth: 7,
                opacity: 0.5,
            }
        },
    });

print('Volumne vs fitted', VAHChart1);

var volDiff = function(feature) {
  var diff = ee.Number(feature.get('vol')).subtract(ee.Number(feature.get('vol_modeled')));
  // Return the feature with the squared difference set to the 'diff' property.
  return feature.set('diff_vol', diff.pow(2));
};

var rmse = ee.Number(
  // Map the difference function over the collection.
  samples.map(volDiff)
  // Reduce to get the mean squared difference.
  .reduceColumns(ee.Reducer.mean(), ['diff_vol'])
  .get('mean')
)
// Compute the square root of the mean square to get RMSE.
.sqrt();

// Print the result.
print('Volume RMSE=', rmse);

/*---------------------------------------------------------------------------------------*/
// Functions

function accumVolume(img,list) {
  // extract out forcing and state variables
  // "past" equivalent to the x(t-1) state of variables
  var past = ee.Image(ee.List(list).get(-1));//.clip(studyArea);
  var pastIt = past.select('Iap');
  var pastPr = past.select('precip');
  var pastAr = past.select('area');
  var pastHt = past.select('height');
  var pastVl = past.select('vol');
  var nowPr = img.select('precip');//.clip(studyArea);
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

// Convert to daily precip
function accumGFS(collection,startDate,nDays) {
  if (nDays>16){
    alert('Max forecast days is 16, only producing forecast for 16 days...');
    nDays = 16;
  }
  var cnt = 1;
  var imgList = [];
  for (var i=0; i<=nDays; i++) {
    var cntMax =(24*(i+1));
    var forecastMeta = [];
    for(cnt;cnt<=cntMax;cnt++){forecastMeta.push(cnt)}
    var dayPrecip = collection.filter(ee.Filter.inList('forecast_hours', forecastMeta));
    imgList.push(dayPrecip.sum().multiply(precipScale)
      .set('system:time_start',startDate.advance(i,'day')));
  }
  return ee.ImageCollection(imgList);
}

function accumChirps(collection, startDate, nDays) {
  // chirps has daily values in it so 
  return ee.ImageCollection(collection.filterDate(startDate, startDate.advance(nDays, 'day')));
}

function timeScale (img){
  return img.multiply(60*60*6)
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
  var prevPrecip = collection.filterDate(s, e); // these lines give you precip of select past days (ie past 7 days)

  var dailyPrev = accumChirps(prevPrecip, s, pastDays);

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


Map.addLayer(pond, {color: 'red'}, 'pond');
Map.centerObject(pond, 13);
