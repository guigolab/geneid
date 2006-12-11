// @(#) $Id: datatype.js,v 1.1 2006-12-11 13:19:07 arnau Exp $

// Run code when the page loads.  From
// http://simon.incutio.com/archive/2004/05/26/addLoadEvent
function addLoadEvent(func) {
  var oldonload = window.onload;
  if (typeof window.onload != 'function') {
    window.onload = func;
  } else {
    window.onload = function() {
      oldonload();
      func();
    }
  }
}

// Set up functions to run when events occur.
function installHandlers() {
  if (!document.getElementById) return;
  var listdatatype = document.getElementById('listdatatype');
  if (listdatatype) {
      // When the user leaves this element, call the server.
      listdatatype.onclick = function() {
          check_datatype(['listdatatype'], ['inputsubmission']);
          return true;          // Continue with default action.
      }
  }
  var fastadatatype = document.getElementById('fastadatatype');
  if (fastadatatype) {
      // When the user leaves this element, call the server.
      fastadatatype.onclick = function() {
          check_datatype(['fastadatatype'], ['inputsubmission']);
          return true;          // Continue with default action.
      }
  }
}

addLoadEvent( installHandlers );
