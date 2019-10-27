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
    mat2 F;
    float mass;
    float volume_0;
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

// simulation parameters
const float dt = 0.1f;
const uint32_t iterations = static_cast<uint32_t>(1.0f / dt);
const vec2 gravity = vec2(0.0f, -0.3f);

// Lamé parameters for stress-strain relationship
const float elastic_lambda = 10.0f;
const float elastic_mu = 20.0f;

uint32_t num_particles = 0;

vector<Particle> ps;
vector<Cell> grid;

WorkerGroup workers;

void clear_grid();
void P2G();
void update_grid();
void G2P();

vector<vec2> spawn_box(int x, int y, int box_x = 16, int box_y = 16)
{
    vector<vec2> temp_positions;
    const float spacing = 0.5f;
    for (float i = - box_x / 2; i < box_x / 2; i += spacing)
    {
        for (float j = - box_y / 2; j < box_y / 2; j += spacing)
        {
            vec2 pos = vec2(x + i, y + j);
            temp_positions.push_back(pos);
        }
    }

    return temp_positions;
}

int main()
{
    vector<vec2> temp_positions = spawn_box(grid_res / 2, grid_res / 2, 16, 16);
    num_particles = temp_positions.size();

    for (int i = 0; i < num_particles; ++i)
    {
        ps.push_back({
            .x = temp_positions[i],
            .v = vec2(rnd() - 0.5f, rnd() - 0.5f + 2.75f) * 0.5f,
            .C = mat2(0.0f),
            .F = mat2(),
            .mass = 1,
        });
    }

    grid.resize(num_cells);

    for (uint32_t i = 0; i < num_cells; i++)
    {
        grid[i].index = i;
    }

    P2G();

    for (auto& p : ps)
    {
        // quadratic interpolation weights
        uvec2 cell_idx = uvec2(p.x);
        vec2 cell_diff = (p.x - vec2(cell_idx)) - vec2(0.5f);

        vec2 weights[3];
        weights[0] = vec2(0.5f) * glm::pow(vec2(0.5f) - cell_diff, vec2(2.0f));
        weights[1] = vec2(0.75f) - glm::pow(cell_diff, vec2(2.0f));
        weights[2] = vec2(0.5f) * glm::pow(vec2(0.5f) + cell_diff, vec2(2.0f));

        float density = 0.0f;
        // iterate over neighbouring 3x3 cells
        for (int gx = 0; gx < 3; ++gx)
        {
            for (int gy = 0; gy < 3; ++gy)
            {
                float weight = weights[gx].x * weights[gy].y;

                // map 2D to 1D index in grid
                int cell_index = ((int)cell_idx.x + (gx - 1)) * grid_res + ((int)cell_idx.y + gy - 1);
                density += grid[cell_index].mass * weight;
            }
        }

        // per-particle volume estimate has now been computed
        float volume = p.mass / density;
        p.volume_0 = volume;
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
            clear_grid();
            P2G();
            update_grid();
            G2P();
        }

        //workers.Run();
        canvas_draw((float*)(&ps[0]), ps.size(), particle_size_bytes);
    };

    emscripten_set_main_loop(loop, 0, 1);

    workers.Terminate();
}

void clear_grid()
{
    for (auto& cell : grid)
    {
        cell.mass = 0;
        cell.v = vec2(0);
    }
}

void P2G()
{
    for (auto& p : ps)
    {
        mat2 F = p.F;

        float J = glm::determinant(F);

        // MPM course, page 46
        float volume = p.volume_0 * J;

        // useful matrices for Neo-Hookean model
        mat2 F_T = glm::transpose(F);
        mat2 F_inv_T = glm::inverse(F_T); // precompute?
        mat2 F_minus_F_inv_T = F - F_inv_T;

        // MPM course equation 48
        mat2 P_term_0 = elastic_mu * (F_minus_F_inv_T);
        mat2 P_term_1 = elastic_lambda * log(J) * F_inv_T;
        mat2 P = P_term_0 + P_term_1;

        // cauchy_stress = (1 / det(F)) * P * F_T
        // equation 38, MPM course
        mat2 stress = (1.0f / J) * (P * F_T);

        // (M_p)^-1 = 4, see APIC paper and MPM course page 42
        // this term is used in MLS-MPM paper eq. 16. with quadratic weights, Mp = (1/4) * (delta_x)^2.
        // in this simulation, delta_x = 1, because i scale the rendering of the domain rather than the domain itself.
        // we multiply by dt as part of the process of fusing the momentum and force update for MLS-MPM
        mat2 eq_16_term_0 = -volume * 4 * stress * dt;

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

                // converting 2D index to 1D
                int cell_index = (int)cell_x.x * grid_res + (int)cell_x.y;
                Cell& cell = grid[cell_index];

                // MPM course, equation 172
                float weighted_mass = weight * p.mass;
                cell.mass += weighted_mass;

                // APIC P2G momentum contribution
                cell.v += weighted_mass * (p.v + Q);

                // fused force/momentum update from MLS-MPM
                // see MLS-MPM paper, equation listed after eqn. 28
                vec2 momentum = (eq_16_term_0 * weight) * cell_dist;
                cell.v += momentum;

                // note: currently "cell.v" refers to MOMENTUM, not velocity!
                // this gets corrected in the UpdateGrid step below.
            }
        }
    }
}

const int wall_min = 3;
const int wall_max = (grid_res - 1) - wall_min;

void update_grid()
{
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
}

void G2P()
{
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

        // NEW: predictive boundary conditions that soften velocities near the domain's edges
        vec2 x_n = p.x + p.v;
        if (x_n.x < wall_min) p.v.x += wall_min - x_n.x;
        if (x_n.x > wall_max) p.v.x += wall_max - x_n.x;
        if (x_n.y < wall_min) p.v.y += wall_min - x_n.y;
        if (x_n.y > wall_max) p.v.y += wall_max - x_n.y;

        // deformation gradient update - MPM course, equation 181
        // Fp' = (I + dt * p.C) * Fp
        mat2 Fp_new = mat2();
        Fp_new += dt * p.C;
        p.F = Fp_new * p.F;
    }
}
