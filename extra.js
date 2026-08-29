// ============================================================
// PERCEPTIONx — pkgdown extra JS
// Sidebars never occupy layout space: they are hidden by CSS and
// pop out as a floating right panel when the toggle button (☰)
// is clicked. Clicking outside the panel closes it again.
// Works with pkgdown >= 2.2 (<main id="main"> + <aside class="col-md-3">).
// ============================================================
(function () {
  "use strict";

  function ready(fn) {
    if (document.readyState !== "loading") {
      fn();
    } else {
      document.addEventListener("DOMContentLoaded", fn);
    }
  }

  ready(function () {
    var aside = document.querySelector("aside");
    if (!aside) return;

    var btn = document.createElement("button");
    btn.id = "toc-toggle-btn";
    btn.type = "button";
    btn.setAttribute("aria-label", "Toggle sidebar");
    btn.setAttribute("title", "Toggle sidebar");
    btn.innerHTML = "&#9776;"; // hamburger

    function setVisible(visible) {
      aside.classList.toggle("toc-visible", visible);
      btn.classList.toggle("active", visible);
    }

    btn.addEventListener("click", function (e) {
      e.stopPropagation();
      setVisible(!aside.classList.contains("toc-visible"));
    });

    // Close the panel when clicking anywhere outside it.
    document.addEventListener("click", function (e) {
      if (aside.classList.contains("toc-visible") &&
          !aside.contains(e.target) && e.target !== btn) {
        setVisible(false);
      }
    });

    // Keep the panel on screen if the page scrolls behind it.
    document.addEventListener("keydown", function (e) {
      if (e.key === "Escape") setVisible(false);
    });

    document.body.appendChild(btn);
  });
})();
