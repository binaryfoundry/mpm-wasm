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

EM_JS(void, canvas_draw, (float* particles, uint32_t size), {
    var ctx = window.ctx;
    var canvas = window.imageCanvas;
    ctx.clearRect(0, 0, canvas.width, canvas.height);
    ctx.beginPath();
    ctx.fillStyle = 'black';
    for (var i = 0; i < size; i++)
    {
        var offset = 8 * i;
        var x = HEAPF32[(particles + offset + 0) >> 2];
        var y = HEAPF32[(particles + offset + 4) >> 2];
        ctx.rect(x, y, 1, 1);
    }
    ctx.fill();
});
