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

const uint32_t screen_width = 640;
const uint32_t screen_height = 640;
const size_t particle_size_bytes = sizeof(Particle);
const uint32_t grid_res = 64;
const uint32_t num_cells = grid_res * grid_res;

uint32_t num_particles = 0;

// simulation parameters

const float dt = 1.0f;
const uint32_t iterations = static_cast<uint32_t>(1.0f / dt);

const vec2 gravity = vec2(0.0f, -0.05f);
vec2 weights[3];

vector<Particle>* ps;
vector<Cell>* grid;

WorkerGroup workers;

int main()
{
    ps = new vector<Particle>();
    grid = new vector<Cell>();

    vector<vec2> temp_positions;
    const float spacing = 1.0f;
    const int box_x = 16, box_y = 16;
    const float sx = grid_res / 2.0f, sy = grid_res / 2.0f;
    for (float i = sx - box_x / 2; i < sx + box_x / 2; i += spacing)
    {
        for (float j = sy - box_y / 2; j < sy + box_y / 2; j += spacing)
        {
            vec2 pos = vec2(i, j);
            temp_positions.push_back(pos);
        }
    }

    num_particles = temp_positions.size();

    for (int i = 0; i < num_particles; ++i)
    {
        ps->push_back({
            .x = temp_positions[i],
            .v = vec2(rnd() - 0.5f, rnd() - 0.5f + 2.75f) * 0.5f,
            .mass = 1
        });
    }

    grid->resize(num_cells);

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
        for (auto& p : *ps)
        {
            //p.x = vec2(rnd() * screen_width, rnd() * screen_height);
        }

        workers.Run();
        canvas_draw((float*)(&(*ps)[0]), ps->size(), particle_size_bytes);
    };

    emscripten_set_main_loop(loop, 0, 1);

    workers.Terminate();

    delete ps;
    delete grid;
}
