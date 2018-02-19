<!DOCTYPE html>
<html lang="en-us">
  <head>
    <meta charset="utf-8">
    <meta http-equiv="Content-Type" content="text/html; charset=utf-8">
    <title>MC Simulation</title>
<script src='//code.jquery.com/jquery-3.3.1.min.js' integrity='sha256-FgpCb/KJQlLNfOu91ta32o/NMZxltwRo8QtmkMRdAu8=' crossorigin='anonymous'></script>
<link href="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/css/bootstrap.min.css" rel="stylesheet" integrity="sha384-Gn5384xqQ1aoWXA+058RXPxPg6fy4IWvTNh0E263XmFcJlSAwiGgFAW/dAiS6JXm" crossorigin="anonymous">
<script src="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/js/bootstrap.min.js" integrity="sha384-JZR6Spejh4U02d8jOt6vLEHfe/JQGiRRSQQxSfFWpi1MquVdAyjUar5+76PVCmYl" crossorigin="anonymous"></script>

<meta name="viewport" content="width=device-width, initial-scale=1">
<link rel='stylesheet' href='style.css'>
  </head>
  <body style='padding: 1em;'>
  <div class='container'>
    <figure style="overflow:visible;" id="spinner"><div class="spinner"></div><center style="margin-top:0.5em"><strong>emscripten</strong></center></figure>
    <div class="emscripten" id="status">Downloading...</div>
    <div class="emscripten">
      <progress value="0" max="100" id="progress" hidden=1></progress>  
    </div>
    <h1>Monte Carlo Simulation</h1>
<?php
  $SFMTbits = array_merge(glob("SFMT*.h"), glob("SFMT*.cpp"));
  echo "<p><a href='http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/SFMT/'>SFMT random numbers</a>: ";
  foreach ($SFMTbits as $f) {
    echo "<a href='$f'>$f</a> ";
  }
  echo "</p>";
  $hBits = glob("*.h");
  echo "<p>Headers: ";
  foreach ($hBits as $f) {
    if (substr($f,0,4) == "SFMT") continue;
    echo "<a href='$f'>$f</a> ";
  }
  echo "</p>";
  $cppBits = glob("*.cpp");
  echo "<p>Source: ";
  foreach ($cppBits as $f) {
    if (substr($f,0,4) == "SFMT") continue;
    echo "<a href='$f'>$f</a> ";
  }
  echo "</p>";
  $extra = array_merge(glob("*.idl"), glob("Makefile"));
  echo "<p>Extra: ";
  foreach ($extra as $f) {
    echo "<a href='$f'>$f</a> ";
  }
  echo "</p>";
?>
    <div id='parametersDiv' class='card'>
      <div class='card-header' role='tab' id='headingParameters'><h5 class='mb-0'><a role='button' data-toggle='collapse' href='#collapseParameters' aria-expanded='true' aria-controls='collapseParameters'>Parameters</a> <button type='button' id='btnNewParameters' style='display: none;' class='btn btn-sm btn-secondary'>New parameters</button></h5></div>
      <div id='collapseParameters' class='collapse show' aria-labelledby='headingParameters'>
      <div class='card-body'>
    <div>
      <div class='row'>
        <div class='col-sm-7 col-md-5 col-lg-4 col-xl-3'><label>Number of molecules: <input class='form-control' size='4' id='numAtoms' value='500' style='width: 4rem;'></label></div>
        <div class='col-sm-6 col-md-4 col-lg-3'><label>Temperature: <input class='form-control' id='T' size='3' value='1' style='width: 5em;'></label></div>
        <div class='col-sm-5 col-md-4 col-lg-3 col-xl-2'><label>Density: <input class='form-control' id='density' size='3' value='1' style='width: 5em;'></label></div>
        <div class='col-sm-8 col-md-6 col-lg-4'><label>mu: <input id='mu' class='form-control' size='3' style='width: 5rem;' disabled></label> <label style='margin-left: 1rem;'><input type='checkbox' id='grandCB'> Grand Canonical</label></div>
      </div>
    <label>Potential: <select class='form-control' id='potType'>
          <option value='0'>Soft sphere</option>
          <option value='1' selected>Lennard-Jones</option>
          <option value='2'>WCA</option>
          <option value='3'>Hard sphere</option>
          <option value='4'>Custom</option>
      </select></label> <span id='ssPowSpan' style='margin-left: 1rem; display: none'><label>Power: <input class='form-control' id='ssPow' size='3' value='12' style='width: 4rem;'></label></span><br>
      <div id='uCustomDiv' style='padding-left: 1rem; display: none'><label>u(r2) = <input class='form-control' id='uCustom' size='40' style='width: 40em;' value='4 * (Math.pow(r2,-6) - Math.pow(r2,-3))'></label>
        <label>r(du/dr) = <input class='form-control' id='duCustom' size='40' style='width: 40em;' value='4 * (-12*Math.pow(r2,-6) + 6*Math.pow(r2,-3))'></label><br>
        <label>r<sup>2</sup>(d<sup>2</sup>u/dr<sup>2</sup>) = <input class='form-control' id='d2uCustom' size='40' style='width: 40em;' value='4 * (12*13*Math.pow(r2,-6) - 6*7*Math.pow(r2,-3))'></label>
        </div>
      <div class='row'><div class='col-sm-12 col-md-6 col-lg-4'>
    <div id='ttDiv'><label>Truncation type: <select class='form-control' id='truncType'>
        <option value='1'>Simple</option>
        <option value='3'>Shifted</option>
        <option value='4' id='truncType4'>Force-shifted</option>
</select></label></div></div>
<div class='col-sm-12 col-md-4 col-lg-3'>
  <div id='rcDiv'><label>Cutoff: <input class='form-control' id='rc' size='3' value='3' style='width: 5em;'></label></div></div>
</div>
  <label><input type='checkbox' id='doMD' onchange='updateDoMD()'> Actually run MD</label> <label>timestep: <input class='form-control' style='width: 4rem;' id='tStep'></label><br>
  <div id='cellsDiv'><label><input type='checkbox' id='doCells'> Use cell lists</label></div>
  <div id='nbrListDiv'><label><input type='checkbox' id='doNbrList'> Use neighbor lists</label> <label>range: <input class='form-control' style='width: 4rem;' id='nbrRange'></label></div>
    <label>Seed: <input class='form-control' id='seed' style='width: 10em;'></label></div>
    <button type='button' id='btnStart' class='btn btn-sm btn-primary'>Start</button></p>
    <div class='output' id="initOutput"></div>
  </div>
  </div>
  </div>

    <div id='resultsDiv' class='card'>
      <div class='card-header' role='tab' id='headingResults'><h5 class='mb-0'><a role='button' data-toggle='collapse' href='#collapseResults' aria-expanded='true' aria-controls='collapseResults'>Simulation Status</a></h5></div>
      <div id='collapseResults' class='collapse show' aria-labelledby='headingResults'>
        <div class='card-body'>
          <div class='output' id="resultsOutput"></div>
          <p class='statusMD statusMC'>Steps: <span id='stepCount'></span><br>
          Speed: <span id='speed'></span> steps/s</p>
          <p class='statusMD'>Simulation time: <span id='simTime'></span></p>
          <p class='statusMC'>Displacement acceptance: <span id='chi'></span><br></p>
          <p class='statusMC'>Step size: <span id='stepSize'></span></p>
          <p class='statusMC'>Insert/Delete acceptance: <span id='chiID'></span></p>
          <p class='statusMC'><button class='btn btn-sm btn-reset btn-success' id='moveNoTune' disabled>Lock step sizes</button></p>
        </div>
      </div>
    </div>
    <div id='dataDiv' class='card'>
      <div class='card-header' role='tab' id='headingData'><h5 class='mb-0'><a role='button' data-toggle='collapse' href='#collapseData' aria-expanded='true' aria-controls='collapseData'>Data</a></h5></div>
      <div id='collapseData' class='collapse show' aria-labelledby='headingData'>
        <div class='card-body'>
          <p><button type='button' id='btnDataEnergy' class='btn btn-sm btn-info' style='display: none;'>Energy</button>
          <button type='button' id='btnDataTemperature' class='btn btn-sm btn-info' style='display: none;'>Temperature</button>
          <button type='button' id='btnDataPressure' class='btn btn-sm btn-info' style='display: none;'>Pressure</button>
          <button type='button' id='btnDataHMA' class='btn btn-sm btn-info' style='display: none;'>HMA</button>
          <button type='button' id='btnDataNA' class='btn btn-sm btn-info' style='display: none;'># of atoms</button>
          <button type='button' id='btnDataDensity' class='btn btn-sm btn-info' style='display: none;'>Density</button></p>
           <div id='dataContent'></div>
        </div>
      </div>
    </div>
  </div>
    <script src='util.js'></script>
    <script src='emscripten.js'></script>
    <script type='text/javascript'>
var box = null, potentialMaster = null, rand = null, move = null, moveID = null, integrator = null, pcHMA = null, meterFull = null, avgFull = null, doMD = false;
var workSteps = 1, totalSteps = 0;
stage = "init";
var running = false, stopRequested = false;
var uCustom = function(r) {return 0;}
var startTime = 0;
function getInputInt(id) {
  var num = document.getElementById(id).value.trim();
  if (!num) return 0;
  num = Number(num);
  //if (num<=0) return 0;
  return num;
}
  function setup() {
    var density = Number(document.getElementById("density").value.trim());
    if (density <= 0) {
      alert("Invalid density "+density);
      return false;
    }
    var potType = Number(document.getElementById("potType").value);
    var truncType = Number(document.getElementById("truncType").value);
    var doLRC = false;
    if (truncType == 2) {
      doLRC = true;
      truncType = 1;
    }
    if (truncType < 0 || truncType > 4 || truncType != Math.round(truncType)) {
      alert("Invalid truncation type " + truncType);
      return false;
    }
    var rc = Number(document.getElementById("rc").value);
    if (truncType>0 && rc <= 0) {
      alert("Invalid rc "+rc);
      document.getElementById("rc").focus();
      return false;
    }
    doMD = document.getElementById("doMD").checked;
    var tStep = doMD ? getInputInt("tStep") : 0;
    var doCells = !doMD && document.getElementById("doCells").checked;
    var doNbrList = doMD && document.getElementById("doNbrList").checked;
    var nbrRange = doNbrList ? getInputInt("nbrRange") : 0;
    if (doMD && !tStep) {
      alert("Must set the timestep");
      document.getElementById("tStep").focus();
      return false;
    }
    if (doNbrList && nbrRange < rc) {
      alert("Neighbor range must be greater than cutoff");
      document.getElementById("nbrRange").focus();
      return false;
    }
    var numAtoms = getInputInt("numAtoms");
    var L = Math.pow(numAtoms/density, 1.0/3.0);
    if ((rc > 0.5*L && !doCells) || rc>L) {
      alert("rc too large, must be less than half the box length ("+L+")");
      document.getElementById("rc").focus();
      return false;
    }

    switch (potType) {
      case 0:
        var ssPow = getInputInt("ssPow");
        potential = new Module.PotentialSS(ssPow, truncType, rc);
        break;
      case 1:
        potential = new Module.PotentialLJ(truncType, rc);
        break;
      case 2:
        potential = new Module.PotentialWCA();
        break;
      case 3:
        potential = new Module.PotentialHS();
        break;
      case 4:
        potential = new PotentialJS();
        eval('potential.u = function(r2) { return '+document.getElementById('uCustom').value+';};');
        eval('potential.du = function(r2) { return '+document.getElementById('duCustom').value+';};');
        eval('potential.d2u = function(r2) { return '+document.getElementById('d2uCustom').value+';};');
        potential.setCutoff(rc);
        potential.setTruncationType(truncType);
        break;
      default:
        throw new Exception("unknown potential type!");
    }
    box = new Module.Box();
    box.setBoxSize(L,L,L);
    box.setNumAtoms(numAtoms);
    box.initCoordinates();
    if (doNbrList) {
      potentialMaster = new Module.PotentialMasterList(potential, box, rc, 2, nbrRange);
    }
    else if (doCells) {
      potentialMaster = new Module.PotentialMasterCell(potential, box, rc, 2);
    }
    else {
      potentialMaster = new Module.PotentialMaster(potential, box);
    }
    if (doCells || doNbrList) potentialMaster.init();
    var seed = getInputInt("seed");
    if (seed == 0) {
      rand = new Module.Random();
      seed = rand.getSeed();
    }
    else {
      rand = new Module.Random(seed);
    }
    document.getElementById("seed").value = seed;
    if (doMD) {
      integrator = new Module.IntegratorMD(potentialMaster, rand, box);
      integrator.setTimeStep(tStep);
      box.enableVelocities();
      if (doNbrList) {
        integrator.setNbrCheckInterval(5);
      }
    }
    else {
      move = new Module.MCMoveDisplacement(box, potentialMaster, rand, 0.2);
      integrator = new Module.IntegratorMC(potentialMaster, rand);
      integrator.addMove(move, 1);
      if (document.getElementById("grandCB").checked) {
        var mu = getInputInt("mu");
        moveID = new Module.MCMoveInsertDelete(box, potentialMaster, rand, mu);
        integrator.addMove(moveID, 1);
      }
    }
    integrator.setTemperature(getInputInt("T"));
    integrator.reset();
    pcHMA = new Module.PotentialCallbackHMA(box, integrator.getTemperature(), 0);

  var fields = ['potType','uCustom','duCustom','d2uCustom','truncType','numAtoms','T','seed','rc','density','mu','grandCB','doCells','doMD','tStep','nbrRange','doNbrList'];
  for (var i in fields) {
    var inp = document.getElementById(fields[i]);
    inp.setAttribute("readonly","true");
    if (inp.type == "checkbox") inp.setAttribute("disabled", "true");
  }
  return true;
}
function startOver() {
  if (running) {
    stopRequested = true;
    running = false;
  }
  steps = 128;
  empty(document.getElementById("stepSizeTable"));
  document.getElementById("stepSizeDiv").style.display = "none";
  $("#collapseStepSize").collapse('show');
  empty(document.getElementById("alphaTable"));
  document.getElementById("alphaDiv").style.display = "none";
  $("#collapseAlpha").collapse('show');
  for (var i=1; i<=nDer; i++) {
    teardownDRow(i);
  }
  teardownDCorTable('Raw');
  teardownDCorTable('Full');
  document.getElementById("productionDiv").style.display = "none";

  var fields = ['potType','uCustom','sigmaHSRef','truncType','NOP','T','seed','nDer','rc'];
  for (var i in fields) {
    document.getElementById(fields[i]).removeAttribute("readonly");
  }
  document.getElementById('setButton').style.display = "";
  var btn = document.getElementById('stepSizeButton');
  btn.style.display = "";
  empty(btn);
  makeText("Find step size", btn);
  btn = document.getElementById('alphaButton');
  btn.style.display = "";
  empty(btn);
  makeText("Find alpha");
  var stages = ['init','stepSize','alpha','production'];
  for (var i=0; i<stages.length; i++) {
    empty(document.getElementById(stages[i]+"Output"));
  }
  stage = "init";
}
document.getElementById('btnStart').addEventListener('click', function(){
  if (running) {
    startTime -= Date.now();
    stopRequested = true;
    return;
  }
  if (box==null) {
    var rv = setup();
    if (!rv) return;
    var doGC = document.getElementById("grandCB").checked;
    var meters = ["Energy", "Pressure", "HMA"];
    if (doGC) {
      meters.push("NA");
      meters.push("Density");
    }
    var doMD = document.getElementById("doMD").checked;
    if (doMD) {
      meters.push("Temperature");
    }
    for (var i=0; i<meters.length; i++) {
      document.getElementById("btnData"+meters[i]).style.display = "";
    }
    document.getElementById("moveNoTune").removeAttribute("disabled");
  }
  startTime += Date.now();
  var mcFunc = function() {integrator.doSteps(workSteps); totalSteps += workSteps;}
  btnRunner("results", "btnStart", mcFunc, updateResults)();
});

function updateResults() {
  var numAtoms = box.getNumAtoms();
  document.getElementById("stepCount").textContent = totalSteps;
  var tStep = getInputInt("tStep");
  document.getElementById("simTime").textContent = (tStep*totalSteps).toPrecision(6);
  var speed = totalSteps/(Date.now()-startTime)*1000;
  document.getElementById("speed").textContent = speed.toPrecision(5);
  if (move) {
    var avgChiGC = moveID ? moveID.getAcceptance() : "";
    var newStepSize = move.get_stepSize();
    var avgChi = move.getAcceptance();
    document.getElementById("stepSize").textContent = newStepSize.toPrecision(6);
    document.getElementById("chi").textContent = avgChi.toPrecision(6);
    if (moveID) document.getElementById("chiID").textContent = avgChiGC.toPrecision(6);
  }
  for (var i=0; i<dataStreams.length; i++) {
    var s = dataStreams[i];
    var d = s.avg.getStatistics();
    d = new Module.ArrayUtil(d);
    var idx = ('idx' in s) ? s.idx : 0;
    var fac = ('fac' in s) ? s.fac : 1;
    document.getElementById("data_"+s.name+"_cur").textContent = (fac*d.x2d(idx,0)).toPrecision(8);
    document.getElementById("data_"+s.name+"_avg").textContent = (fac*d.x2d(idx,1)).toPrecision(8);
    document.getElementById("data_"+s.name+"_err").textContent = (fac*d.x2d(idx,2)).toPrecision(4);
    document.getElementById("data_"+s.name+"_cor").textContent = d.x2d(idx,3).toPrecision(4);
  }
  return true;
}
var steps = 128;
function startRunner(btn, cFunc, cb) {
  var runner = function() {
    if (stopRequested) {
      running = false;
      stopRequested = false;
      btn.textContent = "Continue";
      return;
    }
    btn.textContent = "Pause";
    running = true;
    var t1 = Date.now();
    cFunc();
    var t2 = Date.now();
    if (t2 - t1 < 50) {
      workSteps *= 2;
    }
    if (cb) cb();
    window.setTimeout(runner, 0);
  };
  window.setTimeout(runner, 0);
}

function btnRunner(thisStage, btnID, cFunc, cb) {
  return function() {
    stage = thisStage;
    if (stopRequested && running) return;
    if (running) {
      cb(true);
      if (stage=='production') startTime -= Date.now();
      stopRequested = true;
      return;
    }
    if (stage=='production') startTime += Date.now();
    var btn = document.getElementById(btnID);
    btn.textContent = "Pause";
    startRunner(btn, cFunc, cb);
  };
}
function updatePotType() {
  var potType = Number(document.getElementById("potType").value);
  var truncType = Number(document.getElementById("truncType").value);
  document.getElementById("rcDiv").style.display = (truncType==0||potType==2||potType==3) ? "none" : "block";
  document.getElementById("ttDiv").style.display = (potType==2||potType==3) ? "none" : "block";
  document.getElementById("uCustomDiv").style.display = potType==4 ? "block" : "none";
  document.getElementById("ssPowSpan").style.display = potType==0 ? "inline" : "none";
}
function updateTrunc() {
  var truncType = Number(document.getElementById("truncType").value);
  document.getElementById("rcDiv").style.display = truncType==0 ? "none" : "block";
}
function updateGC() {
  var doGC = document.getElementById("grandCB").checked;
  var muInp = document.getElementById("mu");
  if (doGC) {
    muInp.removeAttribute("disabled");
  }
  else {
    muInp.setAttribute("disabled", "true");
  }
}
function updateDoMD() {
  var doMD = document.getElementById("doMD").checked;
  var mcBits = document.getElementsByClassName("statusMC");
  var mdBits = document.getElementsByClassName("statusMD");
  for (var i=0; i<mcBits.length; i++) mcBits[i].style.display = "none";
  for (var i=0; i<mdBits.length; i++) mdBits[i].style.display = "none";
  if (!doMD) for (var i=0; i<mcBits.length; i++) mcBits[i].style.display = "";
  else for (var i=0; i<mdBits.length; i++) mdBits[i].style.display = "";
  if (doMD) {
    document.getElementById("tStep").removeAttribute("disabled");
  }
  else {
    document.getElementById("tStep").setAttribute("disabled", "true");
  }
  document.getElementById("cellsDiv").style.display = doMD ? "none" : "";
  document.getElementById("nbrListDiv").style.display = doMD ? "" : "none";
}
document.getElementById("potType").addEventListener("change", updatePotType);
document.getElementById("truncType").addEventListener("change", updateTrunc);
document.getElementById("grandCB").addEventListener("change", updateGC);
document.getElementById("moveNoTune").addEventListener("click", function() {
  integrator.setTuning(false);
  var btn = document.getElementById("moveNoTune");
  btn.parentNode.removeChild(btn);
});
window.addEventListener("load", updatePotType);
window.addEventListener("load", updateTrunc);
window.addEventListener("load", updateGC);
window.addEventListener("load", updateDoMD);
var dataStreams = [];
function makeDataDiv(name, av) {
  var row = makeElement("DIV", document.getElementById("dataContent"), {id: "data_"+name, className: "row"});
  var label = makeElement("DIV", row, {className: "col-sm-4 col-md-3 col-lg-3 col-xl-2", textContent: name+" "});
  makeElement("BUTTON", label, {className: 'btn btn-sm btn-success btn-reset', onclick: function(){av.reset();}, textContent: "reset"});
  var cur = makeElement("DIV", row, {className: "col-sm-4 col-md-4 col-lg-3 col-xl-2"});
  makeText("cur: ", cur);
  makeElement("SPAN", cur, {id: "data_"+name+"_cur"});
  var avg = makeElement("DIV", row, {className: "col-sm-4 col-md-4 col-lg-3 col-xl-2"});
  makeText("avg: ", avg);
  makeElement("SPAN", avg, {id: "data_"+name+"_avg"});
  var err = makeElement("DIV", row, {className: "col-sm-4 col-md-3 col-lg-2"});
  makeText("error: ", err);
  makeElement("SPAN", err, {id: "data_"+name+"_err"});
  var cor = makeElement("DIV", row, {className: "col-sm-4 col-md-3 col-lg-2"});
  makeText("cor: ", cor);
  makeElement("SPAN", cor, {id: "data_"+name+"_cor"});
  makeElement("HR", document.getElementById("dataContent"));
}
document.getElementById("btnDataEnergy").addEventListener("click", function() {
    document.getElementById("btnDataEnergy").style.display = "none";
    var meter = new Module.MeterPotentialEnergy(integrator);
    var av = new Module.Average(1, 10, 100);
    var pump = new Module.DataPump(meter, doMD?1:10, av);
    integrator.addListener(pump);
    makeDataDiv("energy", av);
    dataStreams.push({name: "energy", avg: av, fac: 1/box.getNumAtoms()});
});
document.getElementById("btnDataTemperature").addEventListener("click", function() {
    document.getElementById("btnDataTemperature").style.display = "none";
    var meter = new Module.MeterKineticEnergy(box);
    meter.setIntegrator(integrator);
    var av = new Module.Average(2, 1, 100);
    var pump = new Module.DataPump(meter, 4, av);
    integrator.addListener(pump);
    makeDataDiv("T", av);
    dataStreams.push({name: "T", avg: av, fac: (2.0/3.0)/box.getNumAtoms(), idx: 0});
    makeDataDiv("total E", av);
    dataStreams.push({name: "total E", avg: av, fac: 1/box.getNumAtoms(), idx: 1});
});
document.getElementById("btnDataPressure").addEventListener("click", function() {
    document.getElementById("btnDataPressure").style.display = "none";
    if (!meterFull) {
      meterFull = new Module.MeterFullCompute(potentialMaster);
      if (doMD) meterFull.setDoCompute(false);
    }
    var pcp = new Module.PotentialCallbackPressure(box, integrator.getTemperature());
    if (doMD) integrator.addPotentialCallback(pcp, 4);
    var nData0 = meterFull.getNumData();
    meterFull.addCallback(pcp);
    if (!avgFull) {
      avgFull = new Module.Average(1, 1, 100);
      var pump = new Module.DataPump(meterFull, doMD ? 4 : 4*box.getNumAtoms(), avgFull);
      integrator.addListener(pump);
    }
    else {
      avgFull.setNumData(nData0+1);
    }
    makeDataDiv("pressure", avgFull);
    dataStreams.push({name: "pressure", avg: avgFull, idx: nData0});
});
document.getElementById("btnDataHMA").addEventListener("click", function() {
    document.getElementById("btnDataHMA").style.display = "none";
    if (!meterFull) {
      meterFull = new Module.MeterFullCompute(potentialMaster);
      if (doMD) meterFull.setDoCompute(false);
    }
    var nData0 = meterFull.getNumData();
    meterFull.addCallback(pcHMA);
    if (doMD) integrator.addPotentialCallback(pcHMA, 4);
    if (!avgFull) {
      avgFull = new Module.Average(2, 1, 100);
      var pump = new Module.DataPump(meterFull, doMD ? 4 : 4*box.getNumAtoms(), avgFull);
      integrator.addListener(pump);
    }
    else {
      avgFull.setNumData(nData0+2);
    }
    makeDataDiv("HMA U", avgFull);
    dataStreams.push({name: "HMA U", avg: avgFull, idx: nData0, fac: 1/box.getNumAtoms()});
    makeDataDiv("HMA P", avgFull);
    dataStreams.push({name: "HMA P", avg: avgFull, idx: nData0+1});
});
document.getElementById("btnDataNA").addEventListener("click", function() {
    document.getElementById("btnDataNA").style.display = "none";
    var meter = new Module.MeterNumAtoms(box);
    var av = new Module.Average(1, 1000, 100);
    var pump = new Module.DataPump(meter, 10, av);
    integrator.addListener(pump);
    makeDataDiv("# of atoms", av);
    dataStreams.push({name: "# of atoms", avg: av});
});
document.getElementById("btnDataDensity").addEventListener("click", function() {
    if (integrator==null) return;
    document.getElementById("btnDataDensity").style.display = "none";
    var meter = new Module.MeterDensity(box);
    var av = new Module.Average(1, 1000, 100);
    var pump = new Module.DataPump(meter, 10, av);
    integrator.addListener(pump);
    makeDataDiv("density", av);
    dataStreams.push({name: "density", avg: av});
});
document.getElementById("grandCB").removeAttribute("readonly");
document.getElementById("grandCB").removeAttribute("disabled");
document.getElementById("doCells").removeAttribute("readonly");
document.getElementById("doCells").removeAttribute("disabled");
document.getElementById("doMD").removeAttribute("readonly");
document.getElementById("doMD").removeAttribute("disabled");
document.getElementById("doNbrList").removeAttribute("readonly");
document.getElementById("doNbrList").removeAttribute("disabled");
    </script>
    <script async type="text/javascript" src="mc.js"></script>
  </body>
</html>
