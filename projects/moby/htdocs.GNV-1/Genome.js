//
// Genome.js
//
// $Id: Genome.js,v 1.1 2006-12-11 13:19:07 arnau Exp $
//

// ### JavaScript Error Checking ### 
onerror = myOnError;
msgArray = new Array();
urlArray = new Array();
lnoArray = new Array();
function myOnError(msg, url, lno) {
  msgArray[msgArray.length] = msg;
  urlArray[urlArray.length] = url;
  lnoArray[lnoArray.length] = lno;
  return true;
}
function displayErrors() {
  if (msgArray.length == 0) {
    // document.write('<tr><td class="path">This page is being displayed on <B>' + navName + '</B> ' + brNum + ' [NS4==' + NS4 + '] </td></tr>');
    return true;
  }
  for (var i=0; i < msgArray.length; i++) {
    document.writeln('<tr><td class="path"><TABLE width=400 border=1 cellpadding=1 cellspacing=0>');
    document.writeln('<TR><TH class="path" width=50 align=right>Error in file:&nbsp;&nbsp;</TH><TD class="path"> ' + urlArray[i] + '</TD></TR>');
    document.writeln('<TR><TH class="path" width=50 align=right>Line number:&nbsp;&nbsp;</TH><TD class="path"> ' + lnoArray[i] + '</TD></TR>');
    document.writeln('<TR><TH class="path" width=50 align=right>Message:&nbsp;&nbsp;</TH><TD class="path"> ' + msgArray[i] + '</TD></TR>');
    document.writeln('</TABLE><BR></td></tr>');
  }
  return true;
}
function reportErrors() {
  if (msgArray.length == 0) {
    // document.write('<tr><td class="path">This page is being displayed on <B>' + navName + '</B> ' + brNum + ' [NS4==' + NS4 + '] </td></tr>');
    return true;
  }
  for (var i=0; i < msgArray.length; i++) {
    document.writeln('<table class="trailer" border=0 cellpadding=0 cellspacing=0 width="100%">');
    document.writeln('<tr><td class="path">');
    document.writeln('<TABLE width=400 border=1 cellpadding=1 cellspacing=0>');
    document.writeln('<TR><TH class="path" width=50 align=right>Error in file:&nbsp;&nbsp;</TH><TD class="path"> ' + urlArray[i] + '</TD></TR>');
    document.writeln('<TR><TH class="path" width=50 align=right>Line number:&nbsp;&nbsp;</TH><TD class="path"> ' + lnoArray[i] + '</TD></TR>');
    document.writeln('<TR><TH class="path" width=50 align=right>Message:&nbsp;&nbsp;</TH><TD class="path"> ' + msgArray[i] + '</TD></TR>');
    document.writeln('</TABLE></td></tr>');
    document.writeln('</table>');
  }
  return true;
}

// This script and many more are available free online at
// The JavaScript Source!! http://javascript.internet.com
// Jeff Lance (jflance@aol.com)
//-- Begin borrowed
var navName = navigator.appName ;
var brVer = navigator.userAgent;
var brNum;
var reg = new RegExp('/');
function verNumIE() {
  var brVerId = brVer.indexOf('MSIE');
  brNum = brVer.substr(brVerId,8);
}
function verNumOt() {
  var brVerId = brVer.search(reg);
  brNum = brVer.substring(brVerId+1);
}
if (navName == 'Microsoft Internet Explorer') { verNumIE() ; }
else { verNumOt() ; }
//-- End borrowed

// ### Window Size related###
// window.location.href = 'index.html';
if (window != top) { top.location.href = location.href; }
//
var NS4 = 0; 
var bVer = parseFloat(brNum);
if ((bVer < 4.1) && (navName == 'Netscape')) NS4 = 1; // 4.04
function init() {
  if (NS4) { setTimeout("window.onresize = redo", 1000); }
  else { window.location.reload(); }
}
function redo() {
  window.location.reload(); 
}
// window.onload = "init();";
window.onreload = "init();";
window.onresize = "init();";
function onResize() {
  window.onreload = "init();";
  window.onresize = "init();";
}

// ### Misc ###
function displayURL() {
  var sURL = window.location.pathname ; // unescape(window.location.pathname);
  document.write('<B>http://genome.imim.es' + sURL + '</B>');
}
// 
function LastUpdate() {
  var theDate = "";
  theDate = document.lastModified;
  document.write("<B>Last Updated&nbsp;<I>"+theDate+"</I>&nbsp;&copy;&nbsp;Genome BioInformatics Research Laboratory</B>");
}
// 
function LastUpdateSplit() {
  var theDate = "";
  theDate = document.lastModified;
  document.write("<B>Last Updated&nbsp;<I>"+theDate+"</I><BR>&copy;&nbsp;Genome BioInformatics Research Laboratory</B>");
}
// Expiration date: NewIcon("MM/DD/YYYY")
function NewIcon(todate) {
  NewHTML(todate,"&nbsp;<IMG SRC='http://genome.imim.es/g_icons/new_rotating.gif' WIDTH=20 HEIGHT=10 BORDER=0 ALT='NEW'>");
}
function NewHTML(todate,toput) {
  exp = new Date(todate);
  cur = new Date();
  if (cur.getTime() < exp.getTime()) { document.write("&nbsp;"+toput); }
  else { document.write("&nbsp;"); }
}
