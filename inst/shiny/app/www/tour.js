// First-visit guided tour for PERCEPTION-shiny (vanilla JS, no dependencies).
// Shows a short step-by-step walkthrough the first time the app is opened
// (remembered via localStorage), and can be re-triggered from the "?" icon.
(function () {
  'use strict';

  var SEEN_KEY = 'perception_tour_seen';

  var steps = [
    {
      center: true,
      title: 'Welcome to PERCEPTION-shiny',
      text: 'This tool predicts patient response and resistance to cancer treatment ' +
            'from single-cell transcriptomics. Take 30 seconds to walk through the core workflow.'
    },
    {
      target: '#home-go_demo',
      title: '1. Try it instantly',
      text: 'Click "Load Demo" to load demo data (49 genes \u00d7 400 cells + 2 pretrained models). ' +
            'No uploads needed \u2014 explore every feature right away.'
    },
    {
      target: 'a[data-value="data"]',
      title: '2. Data',
      text: 'Upload your own expression matrix, clinical response, and clone annotations, ' +
            'or load the DepMap reference data here.'
    },
    {
      target: 'a[data-value="train"]',
      title: '3. Train',
      text: 'Train drug-response models (elastic net / random forest) on DepMap expression data. ' +
            'Training runs asynchronously in a background process, so the UI stays responsive.'
    },
    {
      target: 'a[data-value="predict"]',
      title: '4. Predict',
      text: 'Use the trained models to predict drug sensitivity for every clone and patient, ' +
            'shown as heatmaps and tables.'
    },
    {
      target: 'a[data-value="visualize"]',
      title: '5. Visualize',
      text: 'ROC curves, response boxplots, clone-viability lollipops, and UMAP spatial views ' +
            'to interpret response and resistance.'
    },
    {
      center: true,
      title: 'You\'re all set!',
      text: 'Click "Load Demo" to start with the demo data, or upload your own. ' +
            'If you want to see this short tutorial again, click the "?" icon in the top navigation bar. Enjoy!'
    }
  ];

  var idx = 0;
  var backdrop = null;
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
  }

  function clearTarget() {
    var t = document.querySelector('.tour-target');
    if (t) t.classList.remove('tour-target');
  }

  function step(i) {
    if (i < 0 || i >= steps.length) return;
    idx = i;
    var s = steps[i];
    clearTarget();

    bubbleTitle.textContent = s.title;
    bubbleText.textContent = s.text;
    document.getElementById('tour-counter').textContent = (i + 1) + ' / ' + steps.length;

    if (s.center) {
      bubble.className = 'active';
      bubble.style.left = '50%';
      bubble.style.top = '38%';
      bubble.style.transform = 'translate(-50%, -50%)';
      backdrop.className = 'active';
    } else {
      var target = document.querySelector(s.target);
      if (!target) { step(i + 1); return; }  // target not on screen: skip quietly
      target.classList.add('tour-target');
      backdrop.className = 'active';

      var r = target.getBoundingClientRect();
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

      if (r.top < 0 || r.bottom > window.innerHeight || r.left < 0 || r.right > window.innerWidth) {
        target.scrollIntoView({ behavior: 'smooth', block: 'center' });
      }
    }

    btnPrev.style.visibility = (idx === 0) ? 'hidden' : 'visible';
    btnNext.textContent = (idx >= steps.length - 1) ? 'Start' : 'Next';
  }

  function end(remember) {
    clearTarget();
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
