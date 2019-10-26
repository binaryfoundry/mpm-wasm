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

EM_JS(void, canvas_draw, (float* particles, uint32_t size, size_t particle_size_bytes), {
    var ctx = window.ctx;
    var canvas = window.imageCanvas;
    ctx.beginPath();
    ctx.fillStyle = '#AED6F1';
    ctx.fillRect(0, 0, canvas.width, canvas.height);
    ctx.fillStyle = 'black';
    var s = (1/64) * 640;
    for (var i = 0; i < size; i++)
    {
        var offset = particle_size_bytes * i;
        var x = HEAPF32[(particles + offset + 0) >> 2];
        var y = HEAPF32[(particles + offset + 4) >> 2];
        ctx.rect((x * s) - 1, 640 - ((y * s) - 1), 3, 3);
    }
    ctx.fill();
});
