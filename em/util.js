
function makeText(str, parentNode) {
  if (!str) return;
  var textNode = document.createTextNode(str);
  if (parentNode) parentNode.appendChild(textNode);
  return textNode;
}

function makeElement(tagName, parentNode, props) {
  var el = document.createElement(tagName);
  if (props) {
    for (p in props) {
      if (p == "disabled" || p=="colspan" || p=="rowspan" || p == "crossorigin" || /^data-/.test(p) || p=="readonly") {
        el.setAttribute(p, props[p]);
      }
      else if (p == "style") {
        for (s in props[p]) {
          el.style[s] = props[p][s];
        }
      }
      else if (p=="onclick" || p=="onchange") {
        if (typeof props[p] == "string") {
          el.setAttribute(p, props[p]);
        }
        else {
          ev = p.replace(/^on/, "");
          el.addEventListener(ev, props[p]);
        }
      }
      else {
        el[p] = props[p];
      }
    }
  }
  if (parentNode) parentNode.appendChild(el);
  return el;
}

function showWarning(warningText) {
  var warning = document.getElementById("warning");
  if (warningText.length==0) {
    warning.style.display = "none";
    return;
  }
  while (warning.firstChild) warning.removeChild(warning.firstChild);
  warning.appendChild(document.createTextNode(warningText));
  warning.style.display = "block";
}

function empty(node) {
  while (node.firstChild) node.removeChild(node.firstChild);
}

var courseDialog = null;
function showCoursesScheduleDialog(courseList) {
  if (!courseDialog) {
    // disable bootstrap buttons
    var bootstrapButton = $.fn.button.noConflict();
    courseDialog = $("#courseDialog").dialog({
      width: 500,
      autoOpen: false,
      buttons: {
        Close: function() {
          courseDialog.dialog( "close" );
        }
      }
    });
    // restore bootstrap buttons
    $.fn.bootstrapBtn = bootstrapButton;
  }
  courseDialog.dialog("open");
  var d = document.getElementById("courseDialog");
  empty(d);
  var cList = makeElement("UL", d);
  for (var i=0; i<courseList.length; i++) {
    var s = makeElement("LI", cList);
    var cData = courseData[courseList[i]];
    var num = cData.cn.replace(/[A-Z]/g, "");
    var l = makeElement("A", s, {style: {color: "blue"}, href: "http://undergrad-catalog.buffalo.edu/courses/index.php?abbr="+cData.subject+"&num="+num, target: "_blank"});
    makeText(cData.subject+" "+cData.cn, l);
    makeText(": "+cData.title, s);
    /*if (dCourses && courseList[i] in dCourses) {
      makeText(" (Diversity)", makeElement("SPAN", s, {style: {color: "#FF4444"}}));
    }*/
    var offList = makeElement("UL", s);
    var typ = makeElement("LI", offList);
    makeText("Typically offered: "+getSemesterList(cData.sMask), typ);
    var sched = makeElement("LI", offList);
    var str = "Scheduled: ";
    if (cData.schedule.length==0) str += "No classes scheduled";
    makeText(str, sched);
    if (cData.schedule.length>0) {
      var sem = getSemesterName(cData.schedule[0]);
      makeText(sem, makeElement("A", sched, {style: {color: "blue"}, target: "_blank", href: "http://www.buffalo.edu/class-schedule?switch=showcourses&semester="+sem.split(" ")[0].toLowerCase()+"&division=UGRD&dept="+cData.subject}));
      for (var j=1; j<cData.schedule.length; j++) {
        makeText(", ", sched);
        var sem = getSemesterName(cData.schedule[j]);
        makeText(sem, makeElement("A", sched, {style: {color: "blue"}, target: "_blank", href: "http://www.buffalo.edu/class-schedule?switch=showcourses&semester="+sem.split(" ")[0].toLowerCase()+"&division=UGRD&dept="+cData.subject}));
      }
    }
  }
}

var semesterNames = {0: "Winter", 1: "Spring", 6: "Summer", 9: "Fall"};
function getSemesterName(term) {
  var y = Math.floor(term/10)+1800;
  return semesterNames[term%10] +" "+ y;
}
function getSemesterList(m) {
  if (m==0) return "Unknown";
  var semesterList = "";
  var isem=0;
  var first = true;
  for (var t in semesterNames) {
    if (m & (1<<isem)) {
      if (!first) semesterList += ", ";
      semesterList += semesterNames[t];
      first = false;
    }
    isem++;
  }
  return semesterList;
}

function validateEmail(email) {
  if (!/^[~{}`\^?=\/*&%$#!'\.A-Za-z_0-9-+]*@[A-Za-z_0-9-\.]*$/.test(email) ||
      /@[^.]*$/.test(email) || /@\./.test(email) || /\.$/.test(email) ||
      /^\./.test(email) || /\.@/.test(email) || /\.\./.test(email)) {
    return false;
  }
  return true;
}

function formatFloatForErr(value, err) {
  if (err==0 || isNaN(err)) return ""+value;
  var r = Math.abs(value/err);
  var rl10 = Math.log10(r);
  var d = rl10+3;
  if (d<3) d=3;
  else if (d>15) d=15;
  return value.toPrecision(d);
}

function getInputInt(id) {
  var num = document.getElementById(id).value.trim();
  if (!num) return 0;
  num = Number(num);
  //if (num<=0) return 0;
  return num;
}
