#include "main.hpp"

#include "canvas.hpp"
#include "worker.hpp"

#include <fstream>
#include <iostream>
#include <vector>

struct Particle
{
    vec2 x;
    vec2 v;
    mat2 C;
    float mass;
    float padding;
};

struct Cell
{
    vec2 v;
    float mass;
    float padding;
};

const uint32_t screen_width = 1280 / 2;
const uint32_t screen_height = 720 / 2;
const uint32_t particle_count = 4096;
const size_t particle_size_bytes = sizeof(Particle);

vector<Particle>* particles;
WorkerGroup workers;

int main()
{
    particles = new vector<Particle>();

    for (uint32_t i = 0; i < particle_count; i++)
    {
        vec2 pos = vec2(rnd() * screen_width, rnd() * screen_height);
        particles->push_back({.x = pos});
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
            p.x = vec2(rnd() * screen_width, rnd() * screen_height);
        }

        workers.Run();
        canvas_draw((float*)(&(*particles)[0]), particles->size(), particle_size_bytes);
    };

    emscripten_set_main_loop(loop, 0, 1);

    workers.Terminate();

    delete particles;
}
