#pragma once

#include <emscripten.h>

EM_JS(void, canvas_setup, (int nx, int ny), {
    var ctx = window.ctx;
    ctx.canvas.width = nx;
    ctx.canvas.height = ny;
    var imageData = ctx.createImageData(nx, ny);
    var imageCanvas = document.createElement('canvas');
    imageCanvas.width = nx;
    imageCanvas.height = ny;
    window.imageData = imageData;
    window.imageCanvas = imageCanvas;
});

EM_JS(void, canvas_draw, (), {
    var ctx = window.ctx;
    var canvas = window.imageCanvas;
    ctx.fillStyle = 'green';
    ctx.fillRect(0, 0, canvas.width, canvas.height);
});
