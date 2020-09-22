uniform sampler2D collided[2]; // single channel
uniform float offset;
uniform float max_idx;
uniform float screenx;
uniform sampler2D barriers;

void main () {
  bool blocked = texture2D(barriers, gl_FragCoord).r > 0;
  float i = gl_FragCoord.x + gl_FragCoord.y * screen.x;
  float idx = mod(i - offset + max_idx, max_idx);
  vec2 coord = vec2(mod(idx, screenx), floor(idx / screenx));
  gl_FragColor = texture2D(blocked ? collided[1] : collided[0], coord);
}