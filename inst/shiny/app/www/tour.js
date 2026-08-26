// First-visit guided tour for PERCEPTION-shiny (vanilla JS, no dependencies).
// Shows a short step-by-step walkthrough the first time the app is opened
// (remembered via localStorage), and can be re-triggered from the "?" icon.
//
// Spotlight approach: a fixed, transparent box (#tour-spotlight) sits over the
// target. It carries BOTH an inset highlight ring AND a huge outset box-shadow
// (0 0 0 9999px) that dims the whole page except the box. Because the box is a
// root-level fixed element, it dims everything regardless of the target's own
// stacking context; and because box-shadow follows border-radius, the hole has
// the SAME rounded corners as the ring. Step targets can request a larger
// square spotlight (e.g. the tiny "?" icon) via the `square` option.
(function () {
  'use strict';

  var SEEN_KEY = 'perception_tour_seen';

  var steps = [
    {
      center: true,
      title: 'Welcome to PERCEPTION-shiny!',
      text: 'This tool predicts patient response and resistance to cancer treatment ' +
            'from single-cell transcriptomics. Take 30 seconds to walk through the core workflow.'
    },
    {
      target: '#home-go_demo',
      title: 'Try it instantly',
      text: 'Click "Load Demo" to load demo data (49 genes × 400 cells + 2 pretrained models). ' +
            'No uploads needed. Explore every feature right away.'
    },
    {
      target: 'a[data-value="data"]',
      title: 'Data',
      text: 'Upload your own expression matrix, clinical response, and clone annotations, ' +
            'or load the DepMap reference data and pretrained models here.'
    },
    {
      target: 'a[data-value="train"]',
      title: 'Train',
      text: 'If the drugs you need are not listed in the pretrained models, try to train a new model here.' +
            'Drug-response models will be trained based on DepMap expression data.'
    },
    {
      target: 'a[data-value="predict"]',
      title: 'Predict',
      text: 'Use the trained models to predict drug sensitivity for every clone and patient, ' +
            'shown as heatmaps and tables.'
    },
    {
      target: 'a[data-value="visualize"]',
      title: 'Visualize',
      text: 'After prediction, you can visualize ROC curves, response boxplots, clone-viability lollipops, and UMAP spatial views ' +
            'to interpret response and resistance.'
    },
    {
      target: '#tour-start',
      box: [48, 56],   // bigger, slightly narrower box centred on the tiny "?" icon
      title: 'See the tutorial again',
      text: 'If you want to see this short tutorial again, click the "?" icon ' +
            'in the top navigation bar.'
    },
    {
      center: true,
      title: 'You\'re all set!',
      text: 'Click "Load Demo" to start with the demo data, or upload your own. Enjoy!'
    }
  ];

  var idx = 0;
  var curStep = null;
  var curTarget = null;
  var backdrop = null;
  var spotlight = null;
  var bubble = null;
  var bubbleTitle = null;
  var bubbleText = null;
  var btnPrev = null;
  var btnNext = null;

  function build() {
    if (document.getElementById('tour-bubble')) return;

    backdrop = document.createElement('div');
    backdrop.id = 'tour-backdrop';
    document.body.appendChild(backdrop);

    spotlight = document.createElement('div');
    spotlight.id = 'tour-spotlight';
    document.body.appendChild(spotlight);

    bubble = document.createElement('div');
    bubble.id = 'tour-bubble';
    bubbleTitle = document.createElement('h4');
    bubbleText = document.createElement('p');

    var counter = document.createElement('span');
    counter.className = 'tour-step-counter';
    counter.id = 'tour-counter';

    btnPrev = document.createElement('button');
    btnPrev.className = 'tour-btn tour-btn-prev';
    btnPrev.textContent = 'Back';
    btnNext = document.createElement('button');
    btnNext.className = 'tour-btn tour-btn-next';
    btnNext.textContent = 'Next';
    var skip = document.createElement('button');
    skip.className = 'tour-btn tour-btn-skip';
    skip.textContent = 'Skip';

    var btns = document.createElement('div');
    btns.className = 'tour-btns';
    btns.appendChild(counter);
    btns.appendChild(btnPrev);
    btns.appendChild(btnNext);
    btns.appendChild(skip);

    bubble.appendChild(bubbleTitle);
    bubble.appendChild(bubbleText);
    bubble.appendChild(btns);
    document.body.appendChild(bubble);

    btnPrev.addEventListener('click', function () { step(idx - 1); });
    btnNext.addEventListener('click', function () {
      if (idx >= steps.length - 1) { end(true); } else { step(idx + 1); }
    });
    skip.addEventListener('click', function () { end(true); });

    // Keep the spotlight glued to the target while scrolling / resizing.
    window.addEventListener('scroll', reposition, { passive: true });
    window.addEventListener('resize', reposition);
  }

  function hideSpotlight() {
    if (spotlight) spotlight.style.display = 'none';
  }

  // Expanded spotlight rectangle for a step: a centred box (w x h) if the step
  // asks for one via `box` / `square`, otherwise the target's own rect.
  function spotlightRect(r, s) {
    if (s.box) {
      var w = s.box[0], h = s.box[1];
      var cx = r.left + r.width / 2;
      var cy = r.top + r.height / 2;
      return {
        left: cx - w / 2, top: cy - h / 2,
        right: cx + w / 2, bottom: cy + h / 2,
        width: w, height: h
      };
    }
    if (s.square) {
      var side = s.square;
      var cx2 = r.left + r.width / 2;
      var cy2 = r.top + r.height / 2;
      return {
        left: cx2 - side / 2, top: cy2 - side / 2,
        right: cx2 + side / 2, bottom: cy2 + side / 2,
        width: side, height: side
      };
    }
    var pad = s.pad || 0;
    return {
      left: r.left - pad, top: r.top - pad,
      right: r.right + pad, bottom: r.bottom + pad,
      width: r.width + 2 * pad, height: r.height + 2 * pad
    };
  }

  function positionSpotlight(rect) {
    spotlight.style.cssText =
      'left:' + rect.left + 'px;top:' + rect.top + 'px;' +
      'width:' + rect.width + 'px;height:' + rect.height + 'px;display:block;';
  }

  function reposition() {
    if (!curTarget || !curStep) return;
    var r = curTarget.getBoundingClientRect();
    if (r.width > 0 && r.height > 0) {
      positionSpotlight(spotlightRect(r, curStep));
    }
  }

  // Place the tooltip bubble near `r` (viewport coordinates), clamped to screen.
  function placeBubble(r) {
    var bw = 340;
    var bh = 170;
    var left = Math.min(Math.max(r.left + r.width / 2 - bw / 2, 12),
                        window.innerWidth - bw - 12);
    var top = r.bottom + 16;
    if (top + bh > window.innerHeight - 12) {
      top = Math.max(r.top - bh - 16, 12);
    }
    bubble.className = 'active';
    bubble.style.left = left + 'px';
    bubble.style.top = top + 'px';
    bubble.style.transform = 'none';
  }

  function step(i) {
    if (i < 0 || i >= steps.length) return;
    idx = i;
    var s = steps[i];
    curStep = s;
    hideSpotlight();
    backdrop.className = '';

    bubbleTitle.textContent = s.title;
    bubbleText.textContent = s.text;
    document.getElementById('tour-counter').textContent = (i + 1) + ' / ' + steps.length;

    if (s.center) {
      curTarget = null;
      bubble.className = 'active';
      bubble.style.left = '50%';
      bubble.style.top = '38%';
      bubble.style.transform = 'translate(-50%, -50%)';
      backdrop.className = 'active';
    } else {
      var target = document.querySelector(s.target);
      if (!target) { step(i + 1); return; }  // target not present: skip quietly
      curTarget = target;
      var r = target.getBoundingClientRect();
      if (r.width === 0 || r.height === 0) {
        // Target is hidden (e.g. collapsed navbar on a narrow viewport) — show
        // the step centered instead of a pointless full-screen dim.
        hideSpotlight();
        bubble.className = 'active';
        bubble.style.left = '50%';
        bubble.style.top = '38%';
        bubble.style.transform = 'translate(-50%, -50%)';
        backdrop.className = 'active';
      } else {
        var sr = spotlightRect(r, s);
        var inView = r.top >= 0 && r.bottom <= window.innerHeight &&
                     r.left >= 0 && r.right <= window.innerWidth;
        if (!inView) {
          target.scrollIntoView({ behavior: 'smooth', block: 'center' });
          // wait for the smooth scroll to settle before placing spotlight
          setTimeout(function () {
            var r2 = target.getBoundingClientRect();
            positionSpotlight(spotlightRect(r2, curStep));
            placeBubble(spotlightRect(r2, curStep));
          }, 450);
        } else {
          positionSpotlight(sr);
          placeBubble(sr);
        }
      }
    }

    btnPrev.style.visibility = (idx === 0) ? 'hidden' : 'visible';
    btnNext.textContent = (idx >= steps.length - 1) ? 'Done' : 'Next';
  }

  function end(remember) {
    curTarget = null;
    curStep = null;
    hideSpotlight();
    if (backdrop) backdrop.className = '';
    if (bubble) bubble.className = '';
    if (remember !== false) {
      try { localStorage.setItem(SEEN_KEY, '1'); } catch (e) { /* private mode */ }
    }
  }

  function start(force) {
    if (!force) {
      try { if (localStorage.getItem(SEEN_KEY)) return; } catch (e) { /* ignore */ }
    }
    build();
    step(0);
  }

  window.addEventListener('load', function () {
    setTimeout(function () { start(false); }, 900);
  });

  // "?" button in the navbar
  document.addEventListener('DOMContentLoaded', function () {
    var btn = document.getElementById('tour-start');
    if (btn) {
      btn.addEventListener('click', function (e) {
        e.preventDefault();
        start(true);
      });
    }
  });

  window.PERCEPTION_TOUR = { start: start };
})();
