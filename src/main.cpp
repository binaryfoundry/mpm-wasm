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
const uint32_t num_particles = 4096;
const size_t particle_size_bytes = sizeof(Particle);

const uint32_t grid_res = 64;
const uint32_t num_cells = grid_res * grid_res;

// simulation parameters

const float dt = 1.0f;
const uint32_t iterations = static_cast<uint32_t>(1.0f / dt);

const vec2 gravity = vec2(0.0f, -0.05f);
vec2 weights[3];

vector<Particle>* particles;
vector<Cell>* cells;

WorkerGroup workers;

int main()
{
    particles = new vector<Particle>();
    cells = new vector<Cell>();

    for (uint32_t i = 0; i < num_particles; i++)
    {
        vec2 pos = vec2(rnd() * screen_width, rnd() * screen_height);
        particles->push_back({
            .x = pos
        });
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
    delete cells;
}
