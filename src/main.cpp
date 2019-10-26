#include <algorithm>
#include <fstream>
#include <iostream>

#include "canvas.hpp"
#include "worker.hpp"

function<void()> update;
void loop() { update(); }

int main()
{
    int threads = std::min((int)thread::hardware_concurrency(), 8);
    cout << "hardware_concurrency: " << threads << std::endl;

    function<void(uint8_t n)> batch_update = [=](uint8_t n) {
    };

    WorkerGroup workers;

    for (int n = 0; n < threads; ++n)
    {
        workers.AddWorker([=] {
            batch_update(n);
        });
    }

    canvas_setup(1280 / 2, 720 / 2);

    update = [&]() {
        workers.Run();
        canvas_draw();
    };

    emscripten_set_main_loop(loop, 0, 1);

    workers.Terminate();
}
