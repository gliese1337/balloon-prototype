export type Lattice = {
  dims: number[];
  vectors: number[][];
  weights: number[];
  opposites: number[];
};

function createVectors(size: number, weights: number[], q: number) {
  const wordsize = q*size;
  const buffer = new Float32Array(wordsize);
  for (let i = 0; i < wordsize; i+=q) {
    buffer.set(weights, i);
  }
  return buffer;
}

function dot(v: number[] | Float32Array, k: number[] | Float32Array, l: number) {
  let d = v[0]*k[0];
  for (let i=1;i<l;i++) d += v[i]*k[i];
  return d;
}

function expandStreamLoop(d: number, q: number, size: number, dims: number[], vectors: number[][], transform: (x: number) => number) {
  const body = [];
  for (let j=1;j<q;j++) {
    const v = vectors[j];
    let offset = v[d-1]*dims[d-2];
    for (let i=d-2;i>0;i--) {
      offset = (v[i]+offset)*dims[i-1];
    }
    offset += v[0];
    if (offset < 0) offset += size;
    const src = transform(j);
    body.push(`streamed[${q}*((i+${offset})%${size})+${j}] = collided[iq${ src < 0 ? ''+src : '+'+src }];`);
  }
  return body.join('\n');
}

function stream(
  d: number, q: number,
  dims: number[], size: number,
  vectors: number[][], opp: number[],
) {
  const body = `
  if (barriers[i]) {
    ${ expandStreamLoop(d, q, size, dims, vectors, x => opp[x]) }
  } else {
    ${ expandStreamLoop(d, q, size, dims, vectors, x => x) }
  }`;

  return new Function('i', 'iq', 'barriers', 'collided', 'streamed', body);
}

function expandPopLoop(d: number, q: number, vectors: number[][]) {
  const body = [];
  for (let j = 1; j < q; j++) {
    const v = vectors[j];
    body.push(`population = streamed[iq+${j}];newrho += population;`);
    for (let k=0;k<d;k++){
      if (v[k] === 0) continue
      body.push(`u${k} ${v[k] < 0 ? '-' : '+'}= population;`);
    }
  }
  return body.join('\n');
}

function expandVectorLoop(d: number, q: number, vectors: number[][], weights: number[]) {
  const body = [];
  for (let j = 1; j < q; j++) {
    const v = vectors[j];
    const terms = [];
    for (let k=0;k<d;k++) {
      if (v[k] === 0) continue;
      terms.push(`${v[k] < 0 ? '-' : '+'}u${k}`);
    }
    
    const direxp = terms.join('');

    body.push(`dir = ${ direxp[0] === '+' ? direxp.substr(1) : direxp }; collided[iq+${j}] = omega * (${weights[j]} * newrho * (usqr + 3 * dir + 4.5 * dir * dir)) + invomega * streamed[iq+${j}]`);
  }

  return body.join('\n');
}

function collide(d: number, q: number, vectors: number[][], weights: number[]) {
  const body = `
    let newrho = streamed[iq];
    let ${ Array.from({ length: d }, (_, i) => `u${i} = 0`).join(', ') };
    let population;
    ${ expandPopLoop(d, q, vectors) }

    if (newrho <= 0) newrho = 0.01;

    rho[i] = newrho;
    ${ Array.from({ length: d }, (_, i) => `u${i} /= newrho;`).join(' ') }

    const tau = 3*viscosity + 0.5; // relaxation timescale
    const omega = 1 / tau;
    const invomega = 1 - omega;

    const usqr =  1 - 1.5 * (${ Array.from({ length: d }, (_, i) => `u${i}*u${i}`).join('+') });
    streamed[iq] = omega * ${weights[0]} * newrho * usqr + invomega * streamed[iq];
    let dir;
    ${ expandVectorLoop(d, q, vectors, weights) }`;
  
  return new Function(
    'i', 'iq', 'viscosity', 'rho', 'streamed', 'collided',
    body
  );
}

export default class LatticeBoltzmann {
  public rho: Float32Array;       // macroscopic density; cached for rendering
  private streamed: Float32Array; // microscopic densities along each lattice direction
  private collided: Float32Array;
  private size: number;           // total number of lattice points
  private dims: number[];
  private vectors: number[][];
  private weights: number[];
  private d: number;
  private q: number;
  private _collide: Function;
  private _stream: Function;

  constructor(lattice: Lattice) {
    const dims = this.dims = lattice.dims;
    const d = this.d = dims.length;
    const weights = this.weights = lattice.weights;
    const q = this.q = weights.length;
    const size = this.size = dims.reduce((a, n) => a * n, 1);
    const vectors = this.vectors = lattice.vectors;
    this.streamed = createVectors(size, weights, q);
    this.collided = createVectors(size, weights, q);
    this.rho = new Float32Array(size);
    this._collide = collide(d, q, vectors, weights);
    this._stream = stream(d, q, dims, size, vectors, lattice.opposites);
  }

  public step(viscosity: number, barriers: boolean[]) {
    const { size, q, _collide, _stream, rho, collided, streamed } = this;

    for (let i=0,iq=0; i<size; i++,iq+=q) {
      _collide(i, iq, viscosity, rho, streamed, collided);
    }
    
    for (let i=0,iq=0; i<size; i++,iq+=q) {
      _stream(i, iq, barriers, collided, streamed);
    }
  }

  // Set all densities in a cell to their equilibrium values for a given velocity and density:
  public setEquilibrium(coords: number[], u: number[], rho: number) {
    const { dims, vectors, weights, d, q, streamed } = this;
  
    let i = coords[d-1];
    for (let j=d-2; j>=0; j--){
      i = i*dims[j] + coords[j];
    }

    this.rho[i] = rho;

    const iq = i*q;
    const usqr =  1 - 1.5 * dot(u, u, d);
    for (let j = 0; j < q; j++) {
      const dir = dot(vectors[j], u, d);
      streamed[iq+j] = weights[j] * rho * (usqr + 3 * dir + 4.5 * dir * dir);
    }
  }
}