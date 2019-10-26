#include "main.hpp"

#include "canvas.hpp"
#include "worker.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>
#include <random>

using std::default_random_engine;
using std::uniform_real_distribution;

default_random_engine generator;
uniform_real_distribution<double> distribution(0.0, 1.0);
inline float rnd() { return static_cast<float>(distribution(generator)); }

function<void()> update;
void loop() { update(); }

struct Particle
{
    vec2 position;
};

const uint32_t screen_width = 1280 / 2;
const uint32_t screen_height = 720 / 2;
const uint32_t particle_count = 4096;

vector<Particle>* particles;
WorkerGroup workers;

int main()
{
    particles = new vector<Particle>();

    for (uint32_t i = 0; i < particle_count; i++)
    {
        vec2 pos = vec2(rnd() * screen_width, rnd() * screen_height);
        particles->push_back({.position = pos});
    }

    int threads = std::min((int)thread::hardware_concurrency(), 8);
    cout << "hardware_concurrency: " << threads << std::endl;

    function<void(uint8_t n)> batch_update = [=](uint8_t n) {
    };

    for (int n = 0; n < threads; ++n)
    {
        workers.AddWorker([=] {
            batch_update(n);
        });
    }

    canvas_setup(screen_width, screen_height);

    update = [&]() {
        for (auto& p : *particles)
        {
            p.position = vec2(rnd() * screen_width, rnd() * screen_height);
        }

        workers.Run();
        canvas_draw((float*)(&(*particles)[0]), particles->size());
    };

    emscripten_set_main_loop(loop, 0, 1);

    workers.Terminate();

    delete particles;
}
