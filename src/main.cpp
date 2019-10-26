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
    uint32_t index;
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

vector<Particle> ps;
vector<Cell> grid;

WorkerGroup workers;

void simulate();

int main()
{
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
        ps.push_back({
            .x = temp_positions[i],
            .v = vec2(rnd() - 0.5f, rnd() - 0.5f + 2.75f) * 0.5f,
            .mass = 1
        });
    }

    grid.resize(num_cells);

    for (uint32_t i = 0; i < num_cells; i++)
    {
        grid[i].index = i;
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
        for (uint32_t i = 0; i < iterations; ++i)
        {
            simulate();
        }

        //workers.Run();
        canvas_draw((float*)(&ps[0]), ps.size(), particle_size_bytes);
    };

    emscripten_set_main_loop(loop, 0, 1);

    workers.Terminate();
}

void simulate()
{
    for (auto& cell : grid)
    {
        cell.mass = 0;
        cell.v = vec2(0);
    }

    // P2G
    for (auto& p : ps)
    {
        // quadratic interpolation weights
        uvec2 cell_idx = uvec2(p.x);
        vec2 cell_diff = (p.x - vec2(cell_idx)) - vec2(0.5f);

        vec2 weights[3];
        weights[0] = vec2(0.5f) * glm::pow(vec2(0.5f) - cell_diff, vec2(2.0f));
        weights[1] = vec2(0.75f) - glm::pow(cell_diff, vec2(2.0f));
        weights[2] = vec2(0.5f) * glm::pow(vec2(0.5f) + cell_diff, vec2(2.0f));

        // for all surrounding 9 cells
        for (uint32_t gx = 0; gx < 3; ++gx)
        {
            for (uint32_t gy = 0; gy < 3; ++gy)
            {
                float weight = weights[gx].x * weights[gy].y;

                uvec2 cell_x = uvec2(cell_idx.x + gx - 1, cell_idx.y + gy - 1);
                vec2 cell_dist = (vec2(cell_x) - p.x) + vec2(0.5f);
                vec2 Q = p.C * cell_dist;

                // MPM course, equation 172
                float mass_contrib = weight * p.mass;

                // converting 2D index to 1D
                int cell_index = (int)cell_x.x * grid_res + (int)cell_x.y;
                Cell& cell = grid[cell_index];

                // scatter mass to the grid
                cell.mass += mass_contrib;
                cell.v += mass_contrib * (p.v + Q);

                // note: currently "cell.v" refers to MOMENTUM, not velocity!
                // this gets corrected in the UpdateGrid step below.
            }
        }
    }

    // grid velocity update
    for (auto& cell : grid)
    {
        if (cell.mass > 0)
        {
            // convert momentum to velocity, apply gravity
            cell.v /= cell.mass;
            cell.v += dt * gravity;

            // boundary conditions
            int x = cell.index / grid_res;
            int y = cell.index % grid_res;
            if (x < 2 || x > grid_res - 3) cell.v.x = 0;
            if (y < 2 || y > grid_res - 3) cell.v.y = 0;
        }
    }

    // G2P
    for (auto& p : ps)
    {
        // reset particle velocity. we calculate it from scratch each step using the grid
        p.v = vec2(0.0f);

        // quadratic interpolation weights
        uvec2 cell_idx = uvec2(p.x);
        vec2 cell_diff = (p.x - vec2(cell_idx)) - vec2(0.5f);

        vec2 weights[3];
        weights[0] = vec2(0.5f) * glm::pow(vec2(0.5f) - cell_diff, vec2(2.0f));
        weights[1] = vec2(0.75f) - glm::pow(cell_diff, vec2(2.0f));
        weights[2] = vec2(0.5f) * glm::pow(vec2(0.5f) + cell_diff, vec2(2.0f));

        // constructing affine per-particle momentum matrix from APIC / MLS-MPM.
        // see APIC paper (https://web.archive.org/web/20190427165435/https://www.math.ucla.edu/~jteran/papers/JSSTS15.pdf), page 6
        // below equation 11 for clarification. this is calculating C = B * (D^-1) for APIC equation 8,
        // where B is calculated in the inner loop at (D^-1) = 4 is a constant when using quadratic interpolation functions
        mat2 B = mat2(0.0f);
        for (uint32_t gx = 0; gx < 3; ++gx)
        {
            for (uint32_t gy = 0; gy < 3; ++gy)
            {
                float weight = weights[gx].x * weights[gy].y;

                uvec2 cell_x = uvec2(cell_idx.x + gx - 1, cell_idx.y + gy - 1);
                int cell_index = (int)cell_x.x * grid_res + (int)cell_x.y;

                vec2 dist = (vec2(cell_x) - p.x) + vec2(0.5f);
                vec2 weighted_velocity = grid[cell_index].v * weight;

                // APIC paper equation 10, constructing inner term for B
                auto term = mat2(weighted_velocity * dist.x, weighted_velocity * dist.y);

                B += term;

                p.v += weighted_velocity;
            }
        }
        p.C = B * 4.0f;

        // advect particles
        p.x += p.v * dt;

        // safety clamp to ensure particles don't exit simulation domain
        p.x = glm::clamp(p.x, vec2(1), vec2(grid_res - 2));

        // mouse interaction
        //if (mouse_down) {
        //    var dist = p.x - mouse_pos;
        //    if (math.dot(dist, dist) < mouse_radius * mouse_radius) {
        //        float norm_factor = (math.length(dist) / mouse_radius);
        //        norm_factor = math.pow(math.sqrt(norm_factor), 8);
        //        var force = math.normalize(dist) * norm_factor * 0.5f;
        //        p.v += force;
        //    }
        //}
    }
}
