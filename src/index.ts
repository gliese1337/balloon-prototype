import LatticeBoltzmann from './simulation'

// Global variables:
const canvas = document.getElementById('theCanvas') as HTMLCanvasElement;
const context = canvas.getContext('2d') as CanvasRenderingContext2D;
const image = context.createImageData(canvas.width, canvas.height);    // for direct pixel manipulation (faster than fillRect)
for (let i=3; i<image.data.length; i+=4) image.data[i] = 255;          // set all alpha values to opaque

const speedSlider = document.getElementById('speedSlider') as HTMLInputElement;
const viscSlider = document.getElementById('viscSlider') as HTMLInputElement;
const mouseSelect = document.getElementById('mouseSelect') as HTMLSelectElement;
const contrastSlider = document.getElementById('contrastSlider') as HTMLSelectElement;

const xdim = canvas.width;
const ydim = canvas.height;

// boolean array of barrier locations
const barrier: boolean[] = Array.from({ length: xdim*ydim }, () => false);

// Create a simple linear "wall" barrier (intentionally a little offset from center):
const barrierSize = ydim/6;
for (let y=(ydim/2)-barrierSize; y<=(ydim/2)+barrierSize; y++) {
  const x = Math.round(ydim/3);
  barrier[x+y*xdim] = true;
}

// Set up the array of colors for plotting (mimicks matplotlib "jet" colormap):
// (Kludge: Index nColors+1 labels the color used for drawing barriers.)
const nColors = 400; // there are actually nColors+2 colors
const colorList: [number, number, number][] = new Array(nColors+2);
for (let c=0; c<=nColors; c++) {
  if (c < nColors/8) {
    colorList[c] = [
      0,
      0,
      Math.round(255 * (c + nColors/8) / (nColors/4)),
    ];
  } else if (c < 3*nColors/8) {
    colorList[c] = [
      0,
      Math.round(255 * (c - nColors/8) / (nColors/4)),
      255,
    ];
  } else if (c < 5*nColors/8) {
    const r = Math.round(255 * (c - 3*nColors/8) / (nColors/4))
    colorList[c] = [
      r,
      255,
      255 - r,
    ];
  } else if (c < 7*nColors/8) {
    colorList[c] = [
      255,
      Math.round(255 * (7*nColors/8 - c) / (nColors/4)),
      0,
    ];
  } else {
    colorList[c] = [
      Math.round(255 * (9*nColors/8 - c) / (nColors/4)),
      0,
      0,
    ];
  }
}

// barriers are black
colorList[nColors+1] = [0,0,0];

// Simulate function executes a bunch of steps and then schedules another call to itself
function simulate(LB: LatticeBoltzmann) {
  const tlimit = Date.now() + 16;
  setBoundaries(LB);

  // Execute a bunch of time steps:
  let c = 0;
  do {
    LB.step(+viscSlider.value, barrier);
    c++;
  } while(Date.now() < tlimit);
  console.log("Iterations per frame:", c);
  paintCanvas(LB);

  requestAnimationFrame(() => simulate(LB));
}

// Make fluid flow in from the left edge
function setBoundaries(LB: LatticeBoltzmann) {
  const u0 = Number(speedSlider.value);
  for (let y=0; y<ydim; y++) {
    LB.setEquilibrium(0, y, u0, 0, 1);
    LB.setEquilibrium(xdim-1, y, u0, 0, 1);
  }
}

function paintCanvas(LB: LatticeBoltzmann) {
  let cIndex=0;
  const contrast = Math.pow(1.2,Number(contrastSlider.value));
  for (let y=0; y<ydim; y++) {
    for (let x=0; x<xdim; x++) {
      if (barrier[LB.index(x,y)]) {
        cIndex = nColors + 1;  // kludge for barrier color which isn't really part of color map
      } else {
        cIndex = Math.round(nColors * ((LB.rho(x,y)-1)*6*contrast + 0.5));
        if (cIndex < 0) cIndex = 0;
        else if (cIndex > nColors) cIndex = nColors;
      }
      image.data.set(colorList[cIndex], 4*(x + y * xdim));
    }
  }

  context.putImageData(image, 0, 0);
}


let mouseIsDown = false;

canvas.addEventListener('mousedown', (e: MouseEvent) => {
  mouseIsDown = true;
  mousePressDrag(e);
}, false);

canvas.addEventListener('mousemove', (e: MouseEvent) => {
  if (mouseIsDown) {
    mousePressDrag(e);
  }
}, false);

document.body.addEventListener('mouseup', () => {
  mouseIsDown = false;
}, false);  // button release could occur outside canvas

// Handle mouse press or drag:
function mousePressDrag(e: MouseEvent) {
  e.preventDefault();
  const x = Math.floor(e.pageX - canvas.offsetLeft);
  const y = Math.floor(e.pageY - canvas.offsetTop);

  if (mouseSelect.selectedIndex === 0) {
    if ((x > 1) && (x < xdim-2) && (y > 1) && (y < ydim-2)) {
      barrier[x+y*xdim] = true;
    }
  } else {
    barrier[x+y*xdim] = false;
  }
}

// Clear all barriers:
(document.getElementById('clearButton') as HTMLButtonElement).addEventListener('click', () => {
  const max = xdim * ydim;
  for (let i=0; i<max; i++) {
    barrier[i] = false;
  }
});

const LB = new LatticeBoltzmann(xdim, ydim);
/*const u0 = Number(speedSlider.value);
for (let y=0; y<ydim; y++) {
  for (let x = 0; x<xdim; x++) {
    LB.setEquilibrium(x, y, u0, 0, 1);
  }
}*/

requestAnimationFrame(() => simulate(LB));