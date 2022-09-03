#include <iostream>
#include <vector>

#include "CGL/vector2D.h"

#include "mass.h"
#include "rope.h"
#include "spring.h"

namespace CGL {

    Rope::Rope(Vector2D start, Vector2D end, int num_nodes, float node_mass, float k, vector<int> pinned_nodes)
    {
        // TODO (Part 1): Create a rope starting at `start`, ending at `end`, and containing `num_nodes` nodes.

        // Comment-in this part when you implement the constructor
        masses.resize(num_nodes);
        for(int i = 0; i <num_nodes; i++)
        {
            masses[i] = new Mass(start + (end - start) / (num_nodes - 1) * i, node_mass, false);
        }
        for (auto &i : pinned_nodes) {
            masses[i]->pinned = true;
        }
        springs.resize(num_nodes - 1);
        for(int i = 0; i < springs.size(); i++)
        {
            springs[i] = new Spring(masses[i], masses[i + 1], k);
        }
    }

    void Rope::simulateEuler(float delta_t, Vector2D gravity)
    {
        for (auto &s : springs)
        {
            Vector2D relative_pos2d = (s->m2->position - s->m1->position);
            Vector2D force          = s->k * relative_pos2d / relative_pos2d.norm() * (relative_pos2d.norm() - s->rest_length);
            s->m1->forces           += force;
            s->m2->forces           -= force;
            // TODO (Part 2): Use Hooke's law to calculate the force on a node
        }

        for (auto &m : masses)
        {
            if (!m->pinned)
            {
                // TODO (Part 2): Add the force due to gravity, then compute the new velocity and position
                m->forces   += m->mass * gravity;

                // TODO (Part 2): Add global damping

                // Explicit method
                // {
                //     m->position += m->velocity * delta_t;
                //     m->velocity += m->forces / m->mass * delta_t;
                // }

                // Semi-implicit method
                {
                    float k_d   = 0.01;
                    m->forces   += - k_d * m->velocity;
                    m->velocity += m->forces / m->mass * delta_t;
                    m->position += m->velocity * delta_t;
                }
            }

            // Reset all forces on each mass
            m->forces = Vector2D(0, 0);
        }
    }

    void Rope::simulateVerlet(float delta_t, Vector2D gravity)
    {
        for (auto &s : springs)
        {
            // TODO (Part 3): Simulate one timestep of the rope using explicit Verlet （solving constraints)
            Vector2D relative_pos2d = (s->m2->position - s->m1->position);
            Vector2D force          = s->k * relative_pos2d / relative_pos2d.norm() * (relative_pos2d.norm() - s->rest_length);
            s->m1->forces           += force;
            s->m2->forces           -= force;
        }

        for (auto &m : masses)
        {
            if (!m->pinned)
            {
                m->forces   += m->mass * gravity;
                Vector2D temp_position = m->position;
                // TODO (Part 3.1): Set the new position of the rope mass

                // verlet
                //m->position += (m->position - m->last_position) + m->forces / m->mass * delta_t * delta_t;

                // dampling
                m->position += (1 - 0.00005) * (m->position - m->last_position) + m->forces / m->mass * delta_t * delta_t;

                m->last_position = temp_position;
                // TODO (Part 4): Add global Verlet damping
            }
        m->forces = Vector2D(0, 0);
        }
    }
}
