#include "WiggleKitServer.h"
#include "godot_cpp/core/class_db.hpp"
#include "godot_cpp/variant/utility_functions.hpp"

using namespace godot;

WiggleKitServer *WiggleKitServer::singleton = nullptr;

WiggleKitServer *WiggleKitServer::get_singleton() {
	return singleton;
}

WiggleKitServer::WiggleKitServer() {
    // The user should still be able to instantiate own instances of
    // the WiggleKitServer. But the one we instantiate in the extension
    // instantiation will be the singleton.
    if (singleton == nullptr)
    {
        singleton = this;
    }
}

WiggleKitServer::~WiggleKitServer() {
    if (singleton == this)
    {
        singleton = nullptr;
    }
}

void WiggleKitServer::_bind_methods() {
	ClassDB::bind_static_method("WiggleKitServer", D_METHOD("get_singleton"), &WiggleKitServer::get_singleton);
	ClassDB::bind_method(D_METHOD("particle_create", "position", "mass", "gravity", "damping"), &WiggleKitServer::particle_create);
	ClassDB::bind_method(D_METHOD("particle_free", "id"), &WiggleKitServer::particle_free);
	ClassDB::bind_method(D_METHOD("particle_set_position", "id", "position"), &WiggleKitServer::particle_set_position);
	ClassDB::bind_method(D_METHOD("particle_get_position", "id"), &WiggleKitServer::particle_get_position);
	ClassDB::bind_method(D_METHOD("particle_set_mass", "id", "mass"), &WiggleKitServer::particle_set_mass);
	ClassDB::bind_method(D_METHOD("particle_get_mass", "id"), &WiggleKitServer::particle_get_mass);
	ClassDB::bind_method(D_METHOD("particle_set_gravity", "id", "gravity"), &WiggleKitServer::particle_set_gravity);
	ClassDB::bind_method(D_METHOD("particle_get_gravity", "id"), &WiggleKitServer::particle_get_gravity);
	ClassDB::bind_method(D_METHOD("particle_set_damping", "id", "damping"), &WiggleKitServer::particle_set_damping);
	ClassDB::bind_method(D_METHOD("particle_get_damping", "id"), &WiggleKitServer::particle_get_damping);
	ClassDB::bind_method(D_METHOD("particle_set_velocity", "id", "velocity"), &WiggleKitServer::particle_set_velocity);
	ClassDB::bind_method(D_METHOD("particle_get_velocity", "id"), &WiggleKitServer::particle_get_velocity);
	ClassDB::bind_method(D_METHOD("distance_constraint_create", "particle_a", "particle_b", "compliance"), &WiggleKitServer::distance_constraint_create);
	ClassDB::bind_method(D_METHOD("distance_constraint_free", "id"), &WiggleKitServer::distance_constraint_free);
	ClassDB::bind_method(D_METHOD("distance_constraint_set_compliance", "id", "compliance"), &WiggleKitServer::distance_constraint_set_compliance);
	ClassDB::bind_method(D_METHOD("distance_constraint_set_rest_length", "id", "rest_length"), &WiggleKitServer::distance_constraint_set_rest_length);
	ClassDB::bind_method(D_METHOD("distance_constraint_get_rest_length", "id"), &WiggleKitServer::distance_constraint_get_rest_length);
	ClassDB::bind_method(D_METHOD("tetrahedral_constraint_create", "particle_a", "particle_b", "particle_c", "particle_d", "compliance"), &WiggleKitServer::tetrahedral_constraint_create);
	ClassDB::bind_method(D_METHOD("tetrahedral_constraint_free", "id"), &WiggleKitServer::tetrahedral_constraint_free);
	ClassDB::bind_method(D_METHOD("tetrahedral_constraint_set_compliance", "id", "compliance"), &WiggleKitServer::tetrahedral_constraint_set_compliance);
	ClassDB::bind_method(D_METHOD("plane_collider_create", "position", "normal", "bounciness", "friction"), &WiggleKitServer::plane_collider_create);
	ClassDB::bind_method(D_METHOD("plane_collider_free", "id"), &WiggleKitServer::plane_collider_free);
	ClassDB::bind_method(D_METHOD("plane_collider_set_position", "id", "position"), &WiggleKitServer::plane_collider_set_position);
	ClassDB::bind_method(D_METHOD("plane_collider_set_normal", "id", "normal"), &WiggleKitServer::plane_collider_set_normal);
	ClassDB::bind_method(D_METHOD("plane_collider_add_particle", "id", "particle_id"), &WiggleKitServer::plane_collider_add_particle);
	ClassDB::bind_method(D_METHOD("plane_collider_remove_particle", "id", "particle_id"), &WiggleKitServer::plane_collider_remove_particle);
	ClassDB::bind_method(D_METHOD("sphere_collider_create", "position", "radius", "bounciness", "friction", "inside"), &WiggleKitServer::sphere_collider_create, DEFVAL(false));
	ClassDB::bind_method(D_METHOD("sphere_collider_free", "id"), &WiggleKitServer::sphere_collider_free);
	ClassDB::bind_method(D_METHOD("sphere_collider_set_position", "id", "position"), &WiggleKitServer::sphere_collider_set_position);
	ClassDB::bind_method(D_METHOD("sphere_collider_set_radius", "id", "radius"), &WiggleKitServer::sphere_collider_set_radius);
	ClassDB::bind_method(D_METHOD("sphere_collider_add_particle", "id", "particle_id"), &WiggleKitServer::sphere_collider_add_particle);
	ClassDB::bind_method(D_METHOD("sphere_collider_remove_particle", "id", "particle_id"), &WiggleKitServer::sphere_collider_remove_particle);
	ClassDB::bind_method(D_METHOD("box_collider_create", "position", "basis", "size", "bounciness", "friction", "inside"), &WiggleKitServer::box_collider_create, DEFVAL(false));
	ClassDB::bind_method(D_METHOD("box_collider_free", "id"), &WiggleKitServer::box_collider_free);
	ClassDB::bind_method(D_METHOD("box_collider_set_position", "id", "position"), &WiggleKitServer::box_collider_set_position);
	ClassDB::bind_method(D_METHOD("box_collider_set_basis", "id", "basis"), &WiggleKitServer::box_collider_set_basis);
	ClassDB::bind_method(D_METHOD("box_collider_set_size", "id", "size"), &WiggleKitServer::box_collider_set_size);
	ClassDB::bind_method(D_METHOD("box_collider_add_particle", "id", "particle_id"), &WiggleKitServer::box_collider_add_particle);
	ClassDB::bind_method(D_METHOD("box_collider_remove_particle", "id", "particle_id"), &WiggleKitServer::box_collider_remove_particle);
	ClassDB::bind_method(D_METHOD("plane_sdf_create", "position", "normal", "rest_distance", "compliance"), &WiggleKitServer::plane_sdf_create);
	ClassDB::bind_method(D_METHOD("plane_sdf_free", "id"), &WiggleKitServer::plane_sdf_free);
	ClassDB::bind_method(D_METHOD("plane_sdf_set_position", "id", "position"), &WiggleKitServer::plane_sdf_set_position);
	ClassDB::bind_method(D_METHOD("plane_sdf_set_normal", "id", "normal"), &WiggleKitServer::plane_sdf_set_normal);
	ClassDB::bind_method(D_METHOD("plane_sdf_set_rest_distance", "id", "rest_distance"), &WiggleKitServer::plane_sdf_set_rest_distance);
	ClassDB::bind_method(D_METHOD("plane_sdf_set_compliance", "id", "compliance"), &WiggleKitServer::plane_sdf_set_compliance);
	ClassDB::bind_method(D_METHOD("plane_sdf_add_particle", "id", "particle_id"), &WiggleKitServer::plane_sdf_add_particle);
	ClassDB::bind_method(D_METHOD("plane_sdf_remove_particle", "id", "particle_id"), &WiggleKitServer::plane_sdf_remove_particle);
	ClassDB::bind_method(D_METHOD("sphere_sdf_create", "position", "radius", "rest_distance", "compliance"), &WiggleKitServer::sphere_sdf_create);
	ClassDB::bind_method(D_METHOD("sphere_sdf_free", "id"), &WiggleKitServer::sphere_sdf_free);
	ClassDB::bind_method(D_METHOD("sphere_sdf_set_position", "id", "position"), &WiggleKitServer::sphere_sdf_set_position);
	ClassDB::bind_method(D_METHOD("sphere_sdf_set_radius", "id", "radius"), &WiggleKitServer::sphere_sdf_set_radius);
	ClassDB::bind_method(D_METHOD("sphere_sdf_set_rest_distance", "id", "rest_distance"), &WiggleKitServer::sphere_sdf_set_rest_distance);
	ClassDB::bind_method(D_METHOD("sphere_sdf_set_compliance", "id", "compliance"), &WiggleKitServer::sphere_sdf_set_compliance);
	ClassDB::bind_method(D_METHOD("sphere_sdf_add_particle", "id", "particle_id"), &WiggleKitServer::sphere_sdf_add_particle);
	ClassDB::bind_method(D_METHOD("sphere_sdf_remove_particle", "id", "particle_id"), &WiggleKitServer::sphere_sdf_remove_particle);
	ClassDB::bind_method(D_METHOD("line_sdf_create", "position_a", "position_b", "rest_distance", "compliance", "infinite"), &WiggleKitServer::line_sdf_create, DEFVAL(true));
	ClassDB::bind_method(D_METHOD("line_sdf_free", "id"), &WiggleKitServer::line_sdf_free);
	ClassDB::bind_method(D_METHOD("line_sdf_set_position_a", "id", "position_a"), &WiggleKitServer::line_sdf_set_position_a);
	ClassDB::bind_method(D_METHOD("line_sdf_set_position_b", "id", "position_b"), &WiggleKitServer::line_sdf_set_position_b);
	ClassDB::bind_method(D_METHOD("line_sdf_set_infinite", "id", "infinite"), &WiggleKitServer::line_sdf_set_infinite);
	ClassDB::bind_method(D_METHOD("line_sdf_set_rest_distance", "id", "rest_distance"), &WiggleKitServer::line_sdf_set_rest_distance);
	ClassDB::bind_method(D_METHOD("line_sdf_set_compliance", "id", "compliance"), &WiggleKitServer::line_sdf_set_compliance);
	ClassDB::bind_method(D_METHOD("line_sdf_add_particle", "id", "particle_id"), &WiggleKitServer::line_sdf_add_particle);
	ClassDB::bind_method(D_METHOD("line_sdf_remove_particle", "id", "particle_id"), &WiggleKitServer::line_sdf_remove_particle);
	ClassDB::bind_method(D_METHOD("box_sdf_create", "position", "basis", "size", "rest_distance", "compliance"), &WiggleKitServer::box_sdf_create);
	ClassDB::bind_method(D_METHOD("box_sdf_free", "id"), &WiggleKitServer::box_sdf_free);
	ClassDB::bind_method(D_METHOD("box_sdf_set_position", "id", "position"), &WiggleKitServer::box_sdf_set_position);
	ClassDB::bind_method(D_METHOD("box_sdf_set_basis", "id", "basis"), &WiggleKitServer::box_sdf_set_basis);
	ClassDB::bind_method(D_METHOD("box_sdf_set_size", "id", "size"), &WiggleKitServer::box_sdf_set_size);
	ClassDB::bind_method(D_METHOD("box_sdf_set_rest_distance", "id", "rest_distance"), &WiggleKitServer::box_sdf_set_rest_distance);
	ClassDB::bind_method(D_METHOD("box_sdf_set_compliance", "id", "compliance"), &WiggleKitServer::box_sdf_set_compliance);
	ClassDB::bind_method(D_METHOD("box_sdf_add_particle", "id", "particle_id"), &WiggleKitServer::box_sdf_add_particle);
	ClassDB::bind_method(D_METHOD("box_sdf_remove_particle", "id", "particle_id"), &WiggleKitServer::box_sdf_remove_particle);
	ClassDB::bind_method(D_METHOD("set_substeps", "substeps"), &WiggleKitServer::set_substeps);
	ClassDB::bind_method(D_METHOD("get_substeps"), &WiggleKitServer::get_substeps);
	ClassDB::bind_method(D_METHOD("iteration", "delta"), &WiggleKitServer::iteration);
	ClassDB::bind_method(D_METHOD("warmup", "iterations", "particle_ids"), &WiggleKitServer::warmup);
}

uint32_t WiggleKitServer::particle_create(const Vector3 &p_position, float p_mass, const Vector3 &p_gravity, float p_damping) {
	uint32_t index;
	if (free_indices.size() > 0) {
		index = free_indices[free_indices.size() - 1];
		free_indices.remove_at(free_indices.size() - 1);
	} else {
		index = particles.size();
		particles.push_back(Particle());
	}

	Particle &p = particles[index];
	p.position = p_position;
	p.previous_position = p_position;
	p.velocity = Vector3();
	p.inv_mass = p_mass > 0.00001f ? 1.0f / p_mass : 0.0f;
	p.gravity = p_gravity;
	p.damping = p_damping;
	p.active = true;

	return index;
}

void WiggleKitServer::particle_free(uint32_t p_id) {
	ERR_FAIL_UNSIGNED_INDEX(p_id, particles.size());
	if (particles[p_id].active) {
		particles[p_id].active = false;
		free_indices.push_back(p_id);

		// Cleanup distance constraints referencing this particle
		for (DistanceConstraint &c : distance_constraints) {
			if (c.active) {
				if (c.particle_a == p_id) {
					c.particle_a = INVALID_ID;
				    c.active = false;
				}
				if (c.particle_b == p_id) {
					c.particle_b = INVALID_ID;
				    c.active = false;
				}
			}
		}

		// Cleanup tetrahedral constraints referencing this particle
		for (TetrahedralConstraint &c : tetrahedral_constraints) {
			if (c.active) {
				for (int i = 0; i < 4; i++) {
					if (c.particles[i] == p_id) {
						c.particles[i] = INVALID_ID;
						c.active = false;
						break;
					}
				}
			}
		}

		// Cleanup plane colliders referencing this particle
		for (PlaneCollider &pc : plane_colliders) {
			if (pc.active) {
				for (int32_t j = (int32_t)pc.particles.size() - 1; j >= 0; j--) {
					if (pc.particles[j] == p_id) {
						pc.particles.remove_at_unordered(j);
					}
				}
			}
		}

		// Cleanup sphere colliders referencing this particle
		for (SphereCollider &sc : sphere_colliders) {
			if (sc.active) {
				for (int32_t j = (int32_t)sc.particles.size() - 1; j >= 0; j--) {
					if (sc.particles[j] == p_id) {
						sc.particles.remove_at_unordered(j);
					}
				}
			}
		}

		// Cleanup box colliders referencing this particle
		for (BoxCollider &bc : box_colliders) {
			if (bc.active) {
				for (int32_t j = (int32_t)bc.particles.size() - 1; j >= 0; j--) {
					if (bc.particles[j] == p_id) {
						bc.particles.remove_at_unordered(j);
					}
				}
			}
		}

		// Cleanup SDF constraints referencing this particle
		for (PlaneSDF &sdf : plane_sdfs) {
			if (sdf.active) {
				for (int32_t j = (int32_t)sdf.particles.size() - 1; j >= 0; j--) {
					if (sdf.particles[j] == p_id) {
						sdf.particles.remove_at_unordered(j);
					}
				}
			}
		}

		for (SphereSDF &sdf : sphere_sdfs) {
			if (sdf.active) {
				for (int32_t j = (int32_t)sdf.particles.size() - 1; j >= 0; j--) {
					if (sdf.particles[j] == p_id) {
						sdf.particles.remove_at_unordered(j);
					}
				}
			}
		}

		for (LineSDF &sdf : line_sdfs) {
			if (sdf.active) {
				for (int32_t j = (int32_t)sdf.particles.size() - 1; j >= 0; j--) {
					if (sdf.particles[j] == p_id) {
						sdf.particles.remove_at_unordered(j);
					}
				}
			}
		}

		for (BoxSDF &sdf : box_sdfs) {
			if (sdf.active) {
				for (int32_t j = (int32_t)sdf.particles.size() - 1; j >= 0; j--) {
					if (sdf.particles[j] == p_id) {
						sdf.particles.remove_at_unordered(j);
					}
				}
			}
		}
	}
}

void WiggleKitServer::particle_set_position(uint32_t p_id, const Vector3 &p_position) {
	ERR_FAIL_UNSIGNED_INDEX(p_id, particles.size());
	Particle &p = particles[p_id];
	ERR_FAIL_COND(!p.active);
	p.position = p_position;
	p.previous_position = p_position;
}

Vector3 WiggleKitServer::particle_get_position(uint32_t p_id) const {
	ERR_FAIL_UNSIGNED_INDEX_V(p_id, particles.size(), Vector3());
	const Particle &p = particles[p_id];
	ERR_FAIL_COND_V(!p.active, Vector3());
	return p.position;
}

void WiggleKitServer::particle_set_mass(uint32_t p_id, float p_mass) {
	ERR_FAIL_UNSIGNED_INDEX(p_id, particles.size());
	Particle &p = particles[p_id];
	ERR_FAIL_COND(!p.active);
	p.inv_mass = p_mass > 0.00001f ? 1.0f / p_mass : 0.0f;
}

float WiggleKitServer::particle_get_mass(uint32_t p_id) const {
	ERR_FAIL_UNSIGNED_INDEX_V(p_id, particles.size(), 0.0f);
	const Particle &p = particles[p_id];
	ERR_FAIL_COND_V(!p.active, 0.0f);
	return p.inv_mass > 0.0f ? 1.0f / p.inv_mass : 0.0f;
}

void WiggleKitServer::particle_set_gravity(uint32_t p_id, const Vector3 &p_gravity) {
	ERR_FAIL_UNSIGNED_INDEX(p_id, particles.size());
	Particle &p = particles[p_id];
	ERR_FAIL_COND(!p.active);
	p.gravity = p_gravity;
}

Vector3 WiggleKitServer::particle_get_gravity(uint32_t p_id) const {
	ERR_FAIL_UNSIGNED_INDEX_V(p_id, particles.size(), Vector3());
	const Particle &p = particles[p_id];
	ERR_FAIL_COND_V(!p.active, Vector3());
	return p.gravity;
}

void WiggleKitServer::particle_set_damping(uint32_t p_id, float p_damping) {
	ERR_FAIL_UNSIGNED_INDEX(p_id, particles.size());
	Particle &p = particles[p_id];
	ERR_FAIL_COND(!p.active);
	p.damping = p_damping;
}

float WiggleKitServer::particle_get_damping(uint32_t p_id) const {
	ERR_FAIL_UNSIGNED_INDEX_V(p_id, particles.size(), 0.0f);
	const Particle &p = particles[p_id];
	ERR_FAIL_COND_V(!p.active, 0.0f);
	return p.damping;
}

void WiggleKitServer::particle_set_velocity(uint32_t p_id, const Vector3 &p_velocity) {
	ERR_FAIL_UNSIGNED_INDEX(p_id, particles.size());
	Particle &p = particles[p_id];
	ERR_FAIL_COND(!p.active);
	p.velocity = p_velocity;
}

Vector3 WiggleKitServer::particle_get_velocity(uint32_t p_id) const {
	ERR_FAIL_UNSIGNED_INDEX_V(p_id, particles.size(), Vector3());
	const Particle &p = particles[p_id];
	ERR_FAIL_COND_V(!p.active, Vector3());
	return p.velocity;
}


uint32_t WiggleKitServer::distance_constraint_create(uint32_t p_particle_a, uint32_t p_particle_b, float p_compliance) {
	ERR_FAIL_UNSIGNED_INDEX_V(p_particle_a, particles.size(), INVALID_ID);
	ERR_FAIL_UNSIGNED_INDEX_V(p_particle_b, particles.size(), INVALID_ID);
	ERR_FAIL_COND_V(!particles[p_particle_a].active, INVALID_ID);
	ERR_FAIL_COND_V(!particles[p_particle_b].active, INVALID_ID);

	uint32_t index;
	if (free_distance_indices.size() > 0) {
		index = free_distance_indices[free_distance_indices.size() - 1];
		free_distance_indices.remove_at(free_distance_indices.size() - 1);
	} else {
		index = distance_constraints.size();
		distance_constraints.push_back(DistanceConstraint());
	}

	DistanceConstraint &c = distance_constraints[index];
	c.particle_a = p_particle_a;
	c.particle_b = p_particle_b;
	c.rest_length = particles[p_particle_a].position.distance_to(particles[p_particle_b].position);
	c.compliance = p_compliance;
	c.active = true;

	return index;
}

void WiggleKitServer::distance_constraint_free(uint32_t p_id) {
	ERR_FAIL_UNSIGNED_INDEX(p_id, distance_constraints.size());
	if (distance_constraints[p_id].active) {
		distance_constraints[p_id].active = false;
		free_distance_indices.push_back(p_id);
	}
}

void WiggleKitServer::distance_constraint_set_compliance(uint32_t p_id, float p_compliance) {
	ERR_FAIL_UNSIGNED_INDEX(p_id, distance_constraints.size());
	ERR_FAIL_COND(!distance_constraints[p_id].active);
	distance_constraints[p_id].compliance = p_compliance;
}

void WiggleKitServer::distance_constraint_set_rest_length(uint32_t p_id, float p_rest_length) {
	ERR_FAIL_UNSIGNED_INDEX(p_id, distance_constraints.size());
	ERR_FAIL_COND(!distance_constraints[p_id].active);
	distance_constraints[p_id].rest_length = p_rest_length;
}

float WiggleKitServer::distance_constraint_get_rest_length(uint32_t p_id) const {
	ERR_FAIL_UNSIGNED_INDEX_V(p_id, distance_constraints.size(), 0.0f);
	ERR_FAIL_COND_V(!distance_constraints[p_id].active, 0.0f);
	return distance_constraints[p_id].rest_length;
}

uint32_t WiggleKitServer::tetrahedral_constraint_create(uint32_t p_particle_a, uint32_t p_particle_b, uint32_t p_particle_c, uint32_t p_particle_d, float p_compliance) {
	ERR_FAIL_UNSIGNED_INDEX_V(p_particle_a, particles.size(), INVALID_ID);
	ERR_FAIL_COND_V(!particles[p_particle_a].active, INVALID_ID);
	ERR_FAIL_UNSIGNED_INDEX_V(p_particle_b, particles.size(), INVALID_ID);
	ERR_FAIL_COND_V(!particles[p_particle_b].active, INVALID_ID);
	ERR_FAIL_UNSIGNED_INDEX_V(p_particle_c, particles.size(), INVALID_ID);
	ERR_FAIL_COND_V(!particles[p_particle_c].active, INVALID_ID);
	ERR_FAIL_UNSIGNED_INDEX_V(p_particle_d, particles.size(), INVALID_ID);
	ERR_FAIL_COND_V(!particles[p_particle_d].active, INVALID_ID);

	uint32_t index;
	if (free_tetrahedral_indices.size() > 0) {
		index = free_tetrahedral_indices[free_tetrahedral_indices.size() - 1];
		free_tetrahedral_indices.remove_at(free_tetrahedral_indices.size() - 1);
	} else {
		index = tetrahedral_constraints.size();
		tetrahedral_constraints.push_back(TetrahedralConstraint());
	}

	TetrahedralConstraint &c = tetrahedral_constraints[index];
	c.particles[0] = p_particle_a;
	c.particles[1] = p_particle_b;
	c.particles[2] = p_particle_c;
	c.particles[3] = p_particle_d;
	
	const Vector3 &p1 = particles[c.particles[0]].position;
	const Vector3 &p2 = particles[c.particles[1]].position;
	const Vector3 &p3 = particles[c.particles[2]].position;
	const Vector3 &p4 = particles[c.particles[3]].position;

	// Volume = 1/6 * (p2-p1) . ((p3-p1) x (p4-p1))
	c.rest_volume = (1.0f / 6.0f) * (p2 - p1).dot((p3 - p1).cross(p4 - p1));
	c.compliance = p_compliance;
	c.active = true;

	return index;
}

void WiggleKitServer::tetrahedral_constraint_free(uint32_t p_id) {
	ERR_FAIL_UNSIGNED_INDEX(p_id, tetrahedral_constraints.size());
	if (tetrahedral_constraints[p_id].active) {
		tetrahedral_constraints[p_id].active = false;
		free_tetrahedral_indices.push_back(p_id);
	}
}

void WiggleKitServer::tetrahedral_constraint_set_compliance(uint32_t p_id, float p_compliance) {
	ERR_FAIL_UNSIGNED_INDEX(p_id, tetrahedral_constraints.size());
	ERR_FAIL_COND(!tetrahedral_constraints[p_id].active);
	tetrahedral_constraints[p_id].compliance = p_compliance;
}

uint32_t WiggleKitServer::plane_collider_create(const Vector3 &p_position, const Vector3 &p_normal, float p_bounciness, float p_friction) {
	uint32_t index;
	if (free_plane_indices.size() > 0) {
		index = free_plane_indices[free_plane_indices.size() - 1];
		free_plane_indices.remove_at(free_plane_indices.size() - 1);
	} else {
		index = plane_colliders.size();
		plane_colliders.push_back(PlaneCollider());
	}

	PlaneCollider &pc = plane_colliders[index];
	pc.position = p_position;
	pc.normal = p_normal.normalized();
	pc.bounciness = p_bounciness;
	pc.friction = p_friction;
	pc.particles.clear();
	pc.active = true;

	return index;
}

void WiggleKitServer::plane_collider_free(uint32_t p_id) {
	ERR_FAIL_UNSIGNED_INDEX(p_id, plane_colliders.size());
	if (plane_colliders[p_id].active) {
		plane_colliders[p_id].active = false;
		free_plane_indices.push_back(p_id);
	}
}

void WiggleKitServer::plane_collider_set_position(uint32_t p_id, const Vector3 &p_position) {
	ERR_FAIL_UNSIGNED_INDEX(p_id, plane_colliders.size());
	ERR_FAIL_COND(!plane_colliders[p_id].active);
	plane_colliders[p_id].position = p_position;
}

void WiggleKitServer::plane_collider_set_normal(uint32_t p_id, const Vector3 &p_normal) {
	ERR_FAIL_UNSIGNED_INDEX(p_id, plane_colliders.size());
	ERR_FAIL_COND(!plane_colliders[p_id].active);
	plane_colliders[p_id].normal = p_normal.normalized();
}

void WiggleKitServer::plane_collider_add_particle(uint32_t p_id, uint32_t p_particle_id) {
	ERR_FAIL_UNSIGNED_INDEX(p_id, plane_colliders.size());
	ERR_FAIL_COND(!plane_colliders[p_id].active);
	ERR_FAIL_UNSIGNED_INDEX(p_particle_id, particles.size());
	ERR_FAIL_COND(!particles[p_particle_id].active);

	plane_colliders[p_id].particles.push_back(p_particle_id);
}

void WiggleKitServer::plane_collider_remove_particle(uint32_t p_id, uint32_t p_particle_id) {
	ERR_FAIL_UNSIGNED_INDEX(p_id, plane_colliders.size());
	ERR_FAIL_COND(!plane_colliders[p_id].active);
	LocalVector<uint32_t> &pts = plane_colliders[p_id].particles;
	for (int32_t j = (int32_t)pts.size() - 1; j >= 0; j--) {
		if (pts[j] == p_particle_id) {
			pts.remove_at_unordered(j);
			return;
		}
	}
}

uint32_t WiggleKitServer::sphere_collider_create(const Vector3 &p_position, float p_radius, float p_bounciness, float p_friction, bool p_inside) {
	uint32_t index;
	if (free_sphere_indices.size() > 0) {
		index = free_sphere_indices[free_sphere_indices.size() - 1];
		free_sphere_indices.remove_at(free_sphere_indices.size() - 1);
	} else {
		index = sphere_colliders.size();
		sphere_colliders.push_back(SphereCollider());
	}

	SphereCollider &sc = sphere_colliders[index];
	sc.position = p_position;
	sc.radius = p_radius;
	sc.bounciness = p_bounciness;
	sc.friction = p_friction;
	sc.particles.clear();
	sc.inside = p_inside;
	sc.active = true;

	return index;
}

void WiggleKitServer::sphere_collider_free(uint32_t p_id) {
	ERR_FAIL_UNSIGNED_INDEX(p_id, sphere_colliders.size());
	if (sphere_colliders[p_id].active) {
		sphere_colliders[p_id].active = false;
		free_sphere_indices.push_back(p_id);
	}
}

void WiggleKitServer::sphere_collider_set_position(uint32_t p_id, const Vector3 &p_position) {
	ERR_FAIL_UNSIGNED_INDEX(p_id, sphere_colliders.size());
	ERR_FAIL_COND(!sphere_colliders[p_id].active);
	sphere_colliders[p_id].position = p_position;
}

void WiggleKitServer::sphere_collider_set_radius(uint32_t p_id, float p_radius) {
	ERR_FAIL_UNSIGNED_INDEX(p_id, sphere_colliders.size());
	ERR_FAIL_COND(!sphere_colliders[p_id].active);
	sphere_colliders[p_id].radius = p_radius;
}

void WiggleKitServer::sphere_collider_add_particle(uint32_t p_id, uint32_t p_particle_id) {
	ERR_FAIL_UNSIGNED_INDEX(p_id, sphere_colliders.size());
	ERR_FAIL_COND(!sphere_colliders[p_id].active);
	ERR_FAIL_UNSIGNED_INDEX(p_particle_id, particles.size());
	ERR_FAIL_COND(!particles[p_particle_id].active);

	sphere_colliders[p_id].particles.push_back(p_particle_id);
}

void WiggleKitServer::sphere_collider_remove_particle(uint32_t p_id, uint32_t p_particle_id) {
	ERR_FAIL_UNSIGNED_INDEX(p_id, sphere_colliders.size());
	ERR_FAIL_COND(!sphere_colliders[p_id].active);
	LocalVector<uint32_t> &pts = sphere_colliders[p_id].particles;
	for (int32_t j = (int32_t)pts.size() - 1; j >= 0; j--) {
		if (pts[j] == p_particle_id) {
			pts.remove_at_unordered(j);
			return;
		}
	}
}

uint32_t WiggleKitServer::box_collider_create(const Vector3 &p_position, const Basis &p_basis, const Vector3 &p_size, float p_bounciness, float p_friction, bool p_inside) {
	uint32_t index;
	if (free_box_indices.size() > 0) {
		index = free_box_indices[free_box_indices.size() - 1];
		free_box_indices.remove_at(free_box_indices.size() - 1);
	} else {
		index = box_colliders.size();
		box_colliders.push_back(BoxCollider());
	}

	BoxCollider &bc = box_colliders[index];
	bc.position = p_position;
	bc.basis = p_basis;
	bc.size = p_size;
	bc.bounciness = p_bounciness;
	bc.friction = p_friction;
	bc.particles.clear();
	bc.inside = p_inside;
	bc.active = true;

	return index;
}

void WiggleKitServer::box_collider_free(uint32_t p_id) {
	ERR_FAIL_UNSIGNED_INDEX(p_id, box_colliders.size());
	if (box_colliders[p_id].active) {
		box_colliders[p_id].active = false;
		free_box_indices.push_back(p_id);
	}
}

void WiggleKitServer::box_collider_set_position(uint32_t p_id, const Vector3 &p_position) {
	ERR_FAIL_UNSIGNED_INDEX(p_id, box_colliders.size());
	ERR_FAIL_COND(!box_colliders[p_id].active);
	box_colliders[p_id].position = p_position;
}

void WiggleKitServer::box_collider_set_basis(uint32_t p_id, const Basis &p_basis) {
	ERR_FAIL_UNSIGNED_INDEX(p_id, box_colliders.size());
	ERR_FAIL_COND(!box_colliders[p_id].active);
	box_colliders[p_id].basis = p_basis;
}

void WiggleKitServer::box_collider_set_size(uint32_t p_id, const Vector3 &p_size) {
	ERR_FAIL_UNSIGNED_INDEX(p_id, box_colliders.size());
	ERR_FAIL_COND(!box_colliders[p_id].active);
	box_colliders[p_id].size = p_size;
}

void WiggleKitServer::box_collider_add_particle(uint32_t p_id, uint32_t p_particle_id) {
	ERR_FAIL_UNSIGNED_INDEX(p_id, box_colliders.size());
	ERR_FAIL_COND(!box_colliders[p_id].active);
	ERR_FAIL_UNSIGNED_INDEX(p_particle_id, particles.size());
	ERR_FAIL_COND(!particles[p_particle_id].active);

	box_colliders[p_id].particles.push_back(p_particle_id);
}

void WiggleKitServer::box_collider_remove_particle(uint32_t p_id, uint32_t p_particle_id) {
	ERR_FAIL_UNSIGNED_INDEX(p_id, box_colliders.size());
	ERR_FAIL_COND(!box_colliders[p_id].active);
	LocalVector<uint32_t> &pts = box_colliders[p_id].particles;
	for (int32_t j = (int32_t)pts.size() - 1; j >= 0; j--) {
		if (pts[j] == p_particle_id) {
			pts.remove_at_unordered(j);
			return;
		}
	}
}

// --- Plane SDF ---

uint32_t WiggleKitServer::plane_sdf_create(const Vector3 &p_position, const Vector3 &p_normal, float p_rest_distance, float p_compliance) {
	uint32_t index;
	if (free_plane_sdf_indices.size() > 0) {
		index = free_plane_sdf_indices[free_plane_sdf_indices.size() - 1];
		free_plane_sdf_indices.remove_at(free_plane_sdf_indices.size() - 1);
	} else {
		index = plane_sdfs.size();
		plane_sdfs.push_back(PlaneSDF());
	}

	PlaneSDF &sdf = plane_sdfs[index];
	sdf.position = p_position;
	sdf.normal = p_normal.normalized();
	sdf.rest_distance = p_rest_distance;
	sdf.compliance = p_compliance;
	sdf.particles.clear();
	sdf.active = true;

	return index;
}

void WiggleKitServer::plane_sdf_free(uint32_t p_id) {
	ERR_FAIL_UNSIGNED_INDEX(p_id, plane_sdfs.size());
	if (plane_sdfs[p_id].active) {
		plane_sdfs[p_id].active = false;
		plane_sdfs[p_id].particles.clear();
		free_plane_sdf_indices.push_back(p_id);
	}
}

void WiggleKitServer::plane_sdf_set_position(uint32_t p_id, const Vector3 &p_position) {
	ERR_FAIL_UNSIGNED_INDEX(p_id, plane_sdfs.size());
	ERR_FAIL_COND(!plane_sdfs[p_id].active);
	plane_sdfs[p_id].position = p_position;
}

void WiggleKitServer::plane_sdf_set_normal(uint32_t p_id, const Vector3 &p_normal) {
	ERR_FAIL_UNSIGNED_INDEX(p_id, plane_sdfs.size());
	ERR_FAIL_COND(!plane_sdfs[p_id].active);
	plane_sdfs[p_id].normal = p_normal.normalized();
}

void WiggleKitServer::plane_sdf_set_rest_distance(uint32_t p_id, float p_rest_distance) {
	ERR_FAIL_UNSIGNED_INDEX(p_id, plane_sdfs.size());
	ERR_FAIL_COND(!plane_sdfs[p_id].active);
	plane_sdfs[p_id].rest_distance = p_rest_distance;
}

void WiggleKitServer::plane_sdf_set_compliance(uint32_t p_id, float p_compliance) {
	ERR_FAIL_UNSIGNED_INDEX(p_id, plane_sdfs.size());
	ERR_FAIL_COND(!plane_sdfs[p_id].active);
	plane_sdfs[p_id].compliance = p_compliance;
}

void WiggleKitServer::plane_sdf_add_particle(uint32_t p_id, uint32_t p_particle_id) {
	ERR_FAIL_UNSIGNED_INDEX(p_id, plane_sdfs.size());
	ERR_FAIL_COND(!plane_sdfs[p_id].active);
	ERR_FAIL_UNSIGNED_INDEX(p_particle_id, particles.size());
	ERR_FAIL_COND(!particles[p_particle_id].active);

	plane_sdfs[p_id].particles.push_back(p_particle_id);
}

void WiggleKitServer::plane_sdf_remove_particle(uint32_t p_id, uint32_t p_particle_id) {
	ERR_FAIL_UNSIGNED_INDEX(p_id, plane_sdfs.size());
	ERR_FAIL_COND(!plane_sdfs[p_id].active);
	LocalVector<uint32_t> &pts = plane_sdfs[p_id].particles;
	for (int32_t j = (int32_t)pts.size() - 1; j >= 0; j--) {
		if (pts[j] == p_particle_id) {
			pts.remove_at_unordered(j);
			return;
		}
	}
}

// --- Sphere SDF ---

uint32_t WiggleKitServer::sphere_sdf_create(const Vector3 &p_position, float p_radius, float p_rest_distance, float p_compliance) {
	uint32_t index;
	if (free_sphere_sdf_indices.size() > 0) {
		index = free_sphere_sdf_indices[free_sphere_sdf_indices.size() - 1];
		free_sphere_sdf_indices.remove_at(free_sphere_sdf_indices.size() - 1);
	} else {
		index = sphere_sdfs.size();
		sphere_sdfs.push_back(SphereSDF());
	}

	SphereSDF &sdf = sphere_sdfs[index];
	sdf.position = p_position;
	sdf.radius = p_radius;
	sdf.rest_distance = p_rest_distance;
	sdf.compliance = p_compliance;
	sdf.particles.clear();
	sdf.active = true;

	return index;
}

void WiggleKitServer::sphere_sdf_free(uint32_t p_id) {
	ERR_FAIL_UNSIGNED_INDEX(p_id, sphere_sdfs.size());
	if (sphere_sdfs[p_id].active) {
		sphere_sdfs[p_id].active = false;
		sphere_sdfs[p_id].particles.clear();
		free_sphere_sdf_indices.push_back(p_id);
	}
}

void WiggleKitServer::sphere_sdf_set_position(uint32_t p_id, const Vector3 &p_position) {
	ERR_FAIL_UNSIGNED_INDEX(p_id, sphere_sdfs.size());
	ERR_FAIL_COND(!sphere_sdfs[p_id].active);
	sphere_sdfs[p_id].position = p_position;
}

void WiggleKitServer::sphere_sdf_set_radius(uint32_t p_id, float p_radius) {
	ERR_FAIL_UNSIGNED_INDEX(p_id, sphere_sdfs.size());
	ERR_FAIL_COND(!sphere_sdfs[p_id].active);
	sphere_sdfs[p_id].radius = p_radius;
}

void WiggleKitServer::sphere_sdf_set_rest_distance(uint32_t p_id, float p_rest_distance) {
	ERR_FAIL_UNSIGNED_INDEX(p_id, sphere_sdfs.size());
	ERR_FAIL_COND(!sphere_sdfs[p_id].active);
	sphere_sdfs[p_id].rest_distance = p_rest_distance;
}

void WiggleKitServer::sphere_sdf_set_compliance(uint32_t p_id, float p_compliance) {
	ERR_FAIL_UNSIGNED_INDEX(p_id, sphere_sdfs.size());
	ERR_FAIL_COND(!sphere_sdfs[p_id].active);
	sphere_sdfs[p_id].compliance = p_compliance;
}

void WiggleKitServer::sphere_sdf_add_particle(uint32_t p_id, uint32_t p_particle_id) {
	ERR_FAIL_UNSIGNED_INDEX(p_id, sphere_sdfs.size());
	ERR_FAIL_COND(!sphere_sdfs[p_id].active);
	ERR_FAIL_UNSIGNED_INDEX(p_particle_id, particles.size());
	ERR_FAIL_COND(!particles[p_particle_id].active);

	sphere_sdfs[p_id].particles.push_back(p_particle_id);
}

void WiggleKitServer::sphere_sdf_remove_particle(uint32_t p_id, uint32_t p_particle_id) {
	ERR_FAIL_UNSIGNED_INDEX(p_id, sphere_sdfs.size());
	ERR_FAIL_COND(!sphere_sdfs[p_id].active);
	LocalVector<uint32_t> &pts = sphere_sdfs[p_id].particles;
	for (int32_t j = (int32_t)pts.size() - 1; j >= 0; j--) {
		if (pts[j] == p_particle_id) {
			pts.remove_at_unordered(j);
			return;
		}
	}
}

// --- Line SDF ---

uint32_t WiggleKitServer::line_sdf_create(const Vector3 &p_position_a, const Vector3 &p_position_b, float p_rest_distance, float p_compliance, bool p_infinite) {
	uint32_t index;
	if (free_line_sdf_indices.size() > 0) {
		index = free_line_sdf_indices[free_line_sdf_indices.size() - 1];
		free_line_sdf_indices.remove_at(free_line_sdf_indices.size() - 1);
	} else {
		index = line_sdfs.size();
		line_sdfs.push_back(LineSDF());
	}

	LineSDF &sdf = line_sdfs[index];
	sdf.position_a = p_position_a;
	sdf.position_b = p_position_b;
	sdf.rest_distance = p_rest_distance;
	sdf.compliance = p_compliance;
	sdf.particles.clear();
	sdf.infinite = p_infinite;
	sdf.active = true;

	return index;
}

void WiggleKitServer::line_sdf_free(uint32_t p_id) {
	ERR_FAIL_UNSIGNED_INDEX(p_id, line_sdfs.size());
	if (line_sdfs[p_id].active) {
		line_sdfs[p_id].active = false;
		line_sdfs[p_id].particles.clear();
		free_line_sdf_indices.push_back(p_id);
	}
}

void WiggleKitServer::line_sdf_set_position_a(uint32_t p_id, const Vector3 &p_position_a) {
	ERR_FAIL_UNSIGNED_INDEX(p_id, line_sdfs.size());
	ERR_FAIL_COND(!line_sdfs[p_id].active);
	line_sdfs[p_id].position_a = p_position_a;
}

void WiggleKitServer::line_sdf_set_position_b(uint32_t p_id, const Vector3 &p_position_b) {
	ERR_FAIL_UNSIGNED_INDEX(p_id, line_sdfs.size());
	ERR_FAIL_COND(!line_sdfs[p_id].active);
	line_sdfs[p_id].position_b = p_position_b;
}

void WiggleKitServer::line_sdf_set_infinite(uint32_t p_id, bool p_infinite) {
	ERR_FAIL_UNSIGNED_INDEX(p_id, line_sdfs.size());
	ERR_FAIL_COND(!line_sdfs[p_id].active);
	line_sdfs[p_id].infinite = p_infinite;
}

void WiggleKitServer::line_sdf_set_rest_distance(uint32_t p_id, float p_rest_distance) {
	ERR_FAIL_UNSIGNED_INDEX(p_id, line_sdfs.size());
	ERR_FAIL_COND(!line_sdfs[p_id].active);
	line_sdfs[p_id].rest_distance = p_rest_distance;
}

void WiggleKitServer::line_sdf_set_compliance(uint32_t p_id, float p_compliance) {
	ERR_FAIL_UNSIGNED_INDEX(p_id, line_sdfs.size());
	ERR_FAIL_COND(!line_sdfs[p_id].active);
	line_sdfs[p_id].compliance = p_compliance;
}

void WiggleKitServer::line_sdf_add_particle(uint32_t p_id, uint32_t p_particle_id) {
	ERR_FAIL_UNSIGNED_INDEX(p_id, line_sdfs.size());
	ERR_FAIL_COND(!line_sdfs[p_id].active);
	ERR_FAIL_UNSIGNED_INDEX(p_particle_id, particles.size());
	ERR_FAIL_COND(!particles[p_particle_id].active);

	line_sdfs[p_id].particles.push_back(p_particle_id);
}

void WiggleKitServer::line_sdf_remove_particle(uint32_t p_id, uint32_t p_particle_id) {
	ERR_FAIL_UNSIGNED_INDEX(p_id, line_sdfs.size());
	ERR_FAIL_COND(!line_sdfs[p_id].active);
	LocalVector<uint32_t> &pts = line_sdfs[p_id].particles;
	for (int32_t j = (int32_t)pts.size() - 1; j >= 0; j--) {
		if (pts[j] == p_particle_id) {
			pts.remove_at_unordered(j);
			return;
		}
	}
}

// --- Box SDF ---

uint32_t WiggleKitServer::box_sdf_create(const Vector3 &p_position, const Basis &p_basis, const Vector3 &p_size, float p_rest_distance, float p_compliance) {
	uint32_t index;
	if (free_box_sdf_indices.size() > 0) {
		index = free_box_sdf_indices[free_box_sdf_indices.size() - 1];
		free_box_sdf_indices.remove_at(free_box_sdf_indices.size() - 1);
	} else {
		index = box_sdfs.size();
		box_sdfs.push_back(BoxSDF());
	}

	BoxSDF &sdf = box_sdfs[index];
	sdf.position = p_position;
	sdf.basis = p_basis;
	sdf.size = p_size;
	sdf.rest_distance = p_rest_distance;
	sdf.compliance = p_compliance;
	sdf.particles.clear();
	sdf.active = true;

	return index;
}

void WiggleKitServer::box_sdf_free(uint32_t p_id) {
	ERR_FAIL_UNSIGNED_INDEX(p_id, box_sdfs.size());
	if (box_sdfs[p_id].active) {
		box_sdfs[p_id].active = false;
		box_sdfs[p_id].particles.clear();
		free_box_sdf_indices.push_back(p_id);
	}
}

void WiggleKitServer::box_sdf_set_position(uint32_t p_id, const Vector3 &p_position) {
	ERR_FAIL_UNSIGNED_INDEX(p_id, box_sdfs.size());
	ERR_FAIL_COND(!box_sdfs[p_id].active);
	box_sdfs[p_id].position = p_position;
}

void WiggleKitServer::box_sdf_set_basis(uint32_t p_id, const Basis &p_basis) {
	ERR_FAIL_UNSIGNED_INDEX(p_id, box_sdfs.size());
	ERR_FAIL_COND(!box_sdfs[p_id].active);
	box_sdfs[p_id].basis = p_basis;
}

void WiggleKitServer::box_sdf_set_size(uint32_t p_id, const Vector3 &p_size) {
	ERR_FAIL_UNSIGNED_INDEX(p_id, box_sdfs.size());
	ERR_FAIL_COND(!box_sdfs[p_id].active);
	box_sdfs[p_id].size = p_size;
}

void WiggleKitServer::box_sdf_set_rest_distance(uint32_t p_id, float p_rest_distance) {
	ERR_FAIL_UNSIGNED_INDEX(p_id, box_sdfs.size());
	ERR_FAIL_COND(!box_sdfs[p_id].active);
	box_sdfs[p_id].rest_distance = p_rest_distance;
}

void WiggleKitServer::box_sdf_set_compliance(uint32_t p_id, float p_compliance) {
	ERR_FAIL_UNSIGNED_INDEX(p_id, box_sdfs.size());
	ERR_FAIL_COND(!box_sdfs[p_id].active);
	box_sdfs[p_id].compliance = p_compliance;
}

void WiggleKitServer::box_sdf_add_particle(uint32_t p_id, uint32_t p_particle_id) {
	ERR_FAIL_UNSIGNED_INDEX(p_id, box_sdfs.size());
	ERR_FAIL_COND(!box_sdfs[p_id].active);
	ERR_FAIL_UNSIGNED_INDEX(p_particle_id, particles.size());
	ERR_FAIL_COND(!particles[p_particle_id].active);

	box_sdfs[p_id].particles.push_back(p_particle_id);
}

void WiggleKitServer::box_sdf_remove_particle(uint32_t p_id, uint32_t p_particle_id) {
	ERR_FAIL_UNSIGNED_INDEX(p_id, box_sdfs.size());
	ERR_FAIL_COND(!box_sdfs[p_id].active);
	LocalVector<uint32_t> &pts = box_sdfs[p_id].particles;
	for (int32_t j = (int32_t)pts.size() - 1; j >= 0; j--) {
		if (pts[j] == p_particle_id) {
			pts.remove_at_unordered(j);
			return;
		}
	}
}

void WiggleKitServer::set_substeps(uint32_t p_substeps) {
	substeps = p_substeps > 0 ? p_substeps : 1;
}

uint32_t WiggleKitServer::get_substeps() const {
	return substeps;
}

void WiggleKitServer::iteration(float p_delta) {
	if (p_delta <= 0.0f || substeps == 0) {
		return;
	}

	float sub_delta = p_delta / (float)substeps;
	float inv_sub_delta = 1.0f / sub_delta;
	float inv_sub_delta_sq = inv_sub_delta * inv_sub_delta;

	for (uint32_t step = 0; step < substeps; step++) {
		// 1. Predict positions
		for (uint32_t i = 0; i < particles.size(); i++) {
			Particle &p = particles[i];
			if (!p.active || p.inv_mass == 0.0f) {
				continue;
			}

			// Exponential velocity decay for damping
			p.velocity *= std::exp(-p.damping * sub_delta);

			p.previous_position = p.position;
			p.velocity += p.gravity * sub_delta;
			p.position += p.velocity * sub_delta;
		}

		// 2. Constraints (XPBD)
		for (uint32_t i = 0; i < distance_constraints.size(); i++) {
			DistanceConstraint &c = distance_constraints[i];
			if (!c.active || c.particle_a == INVALID_ID || c.particle_b == INVALID_ID) {
				continue;
			}

			Particle &pa = particles[c.particle_a];
			Particle &pb = particles[c.particle_b];

			float w_a = pa.inv_mass;
			float w_b = pb.inv_mass;
			float w_sum = w_a + w_b;

			if (w_sum == 0.0f) {
				continue;
			}

			Vector3 diff = pa.position - pb.position;
			float dist = diff.length();
			if (dist < 0.0001f) {
				continue;
			}

			Vector3 n = diff / dist;
			float C = dist - c.rest_length;
			float alpha = c.compliance * inv_sub_delta_sq;

			// Simplified XPBD for substeps: No lambda tracking across iterations, 
			// because there is only 1 "iteration" per substep.
			float d_lambda = -C / (w_sum + alpha);

			pa.position += n * (d_lambda * w_a);
			pb.position -= n * (d_lambda * w_b);
		}

		for (uint32_t i = 0; i < tetrahedral_constraints.size(); i++) {
			TetrahedralConstraint &c = tetrahedral_constraints[i];
			if (!c.active) continue;

			Particle &p1 = particles[c.particles[0]];
			Particle &p2 = particles[c.particles[1]];
			Particle &p3 = particles[c.particles[2]];
			Particle &p4 = particles[c.particles[3]];

			float w1 = p1.inv_mass;
			float w2 = p2.inv_mass;
			float w3 = p3.inv_mass;
			float w4 = p4.inv_mass;

			if (w1 + w2 + w3 + w4 == 0.0f) continue;

			// Gradients for C = 1/6 * (p2-p1) . ((p3-p1) x (p4-p1)) - V_rest
			Vector3 grad2 = (1.0f / 6.0f) * (p3.position - p1.position).cross(p4.position - p1.position);
			Vector3 grad3 = (1.0f / 6.0f) * (p4.position - p1.position).cross(p2.position - p1.position);
			Vector3 grad4 = (1.0f / 6.0f) * (p2.position - p1.position).cross(p3.position - p1.position);
			Vector3 grad1 = -(grad2 + grad3 + grad4);

			float w_sum = w1 * grad1.length_squared() + w2 * grad2.length_squared() +
						  w3 * grad3.length_squared() + w4 * grad4.length_squared();

			if (w_sum == 0.0f) continue;

			float current_vol = (1.0f / 6.0f) * (p2.position - p1.position).dot((p3.position - p1.position).cross(p4.position - p1.position));
			float C = current_vol - c.rest_volume;
			float alpha = c.compliance * inv_sub_delta_sq;

			float d_lambda = -C / (w_sum + alpha);

			p1.position += d_lambda * w1 * grad1;
			p2.position += d_lambda * w2 * grad2;
			p3.position += d_lambda * w3 * grad3;
			p4.position += d_lambda * w4 * grad4;
		}

		// SDF constraints (XPBD)
		for (PlaneSDF &sdf : plane_sdfs) {
			if (!sdf.active) continue;
			float alpha = sdf.compliance * inv_sub_delta_sq;

			for (uint32_t p_idx : sdf.particles) {
				Particle &p = particles[p_idx];
				if (!p.active || p.inv_mass == 0.0f) continue;

				// Signed distance from plane: C = dot(pos - plane_pos, normal) - rest_distance
				float signed_dist = sdf.normal.dot(p.position - sdf.position);
				float C = signed_dist - sdf.rest_distance;

				// grad C = normal, |grad C|^2 = 1
				float d_lambda = -C / (p.inv_mass + alpha);
				p.position += sdf.normal * (d_lambda * p.inv_mass);
			}
		}

		for (SphereSDF &sdf : sphere_sdfs) {
			if (!sdf.active) continue;
			float alpha = sdf.compliance * inv_sub_delta_sq;

			for (uint32_t p_idx : sdf.particles) {
				Particle &p = particles[p_idx];
				if (!p.active || p.inv_mass == 0.0f) continue;

				Vector3 diff = p.position - sdf.position;
				float dist = diff.length();
				if (dist < 0.0001f) continue;

				// Signed distance from sphere surface: positive = outside
				float signed_dist = dist - sdf.radius;
				float C = signed_dist - sdf.rest_distance;

				Vector3 grad = diff / dist;
				float d_lambda = -C / (p.inv_mass + alpha);
				p.position += grad * (d_lambda * p.inv_mass);
			}
		}

		for (LineSDF &sdf : line_sdfs) {
			if (!sdf.active) continue;
			float alpha = sdf.compliance * inv_sub_delta_sq;

			Vector3 ab = sdf.position_b - sdf.position_a;
			float ab_len_sq = ab.length_squared();

			for (uint32_t p_idx : sdf.particles) {
				Particle &p = particles[p_idx];
				if (!p.active || p.inv_mass == 0.0f) continue;

				// Project onto line to find closest point
				Vector3 ap = p.position - sdf.position_a;
				float t = ab_len_sq > 0.0001f ? ap.dot(ab) / ab_len_sq : 0.0f;

				// Clamp to segment if not infinite
				if (!sdf.infinite) {
					if (t < 0.0f) t = 0.0f;
					else if (t > 1.0f) t = 1.0f;
				}

				Vector3 closest = sdf.position_a + ab * t;
				Vector3 diff = p.position - closest;
				float dist = diff.length();
				if (dist < 0.0001f) continue;

				float C = dist - sdf.rest_distance;

				Vector3 grad = diff / dist;
				float d_lambda = -C / (p.inv_mass + alpha);
				p.position += grad * (d_lambda * p.inv_mass);
			}
		}

		for (BoxSDF &sdf : box_sdfs) {
			if (!sdf.active) continue;
			float alpha = sdf.compliance * inv_sub_delta_sq;

			for (uint32_t p_idx : sdf.particles) {
				Particle &p = particles[p_idx];
				if (!p.active || p.inv_mass == 0.0f) continue;

				// Transform to local box space
				Vector3 local_pos = sdf.basis.xform_inv(p.position - sdf.position);
				Vector3 half_extents = sdf.size * 0.5f;

				// Signed distance to box (negative inside, positive outside)
				Vector3 d = local_pos.abs() - half_extents;
				float outside_dist = Vector3(UtilityFunctions::maxf(d.x, 0.0f), UtilityFunctions::maxf(d.y, 0.0f), UtilityFunctions::maxf(d.z, 0.0f)).length();
				float inside_dist = UtilityFunctions::minf(UtilityFunctions::maxf(d.x, UtilityFunctions::maxf(d.y, d.z)), 0.0f);
				float signed_dist = outside_dist + inside_dist;

				float C = signed_dist - sdf.rest_distance;

				// Compute gradient (normal direction of SDF)
				Vector3 local_grad;
				if (d.x > 0.0f || d.y > 0.0f || d.z > 0.0f) {
					// Outside: gradient points from closest surface point to particle
					Vector3 clamped;
					clamped.x = UtilityFunctions::clampf(local_pos.x, -half_extents.x, half_extents.x);
					clamped.y = UtilityFunctions::clampf(local_pos.y, -half_extents.y, half_extents.y);
					clamped.z = UtilityFunctions::clampf(local_pos.z, -half_extents.z, half_extents.z);
					local_grad = local_pos - clamped;
					float len = local_grad.length();
					if (len < 0.0001f) continue;
					local_grad /= len;
				} else {
					// Inside: gradient points toward nearest face
					int axis = 0;
					float max_d = d.x;
					if (d.y > max_d) { axis = 1; max_d = d.y; }
					if (d.z > max_d) { axis = 2; max_d = d.z; }
					local_grad = Vector3(0, 0, 0);
					local_grad[axis] = local_pos[axis] > 0 ? 1.0f : -1.0f;
				}

				Vector3 grad = sdf.basis.xform(local_grad);
				float d_lambda = -C / (p.inv_mass + alpha);
				p.position += grad * (d_lambda * p.inv_mass);
			}
		}

		// 3. Collision constraints (XPBD-style)
		for (PlaneCollider &pc : plane_colliders) {
			if (!pc.active) {
				continue;
			}

			for (uint32_t p_idx : pc.particles) {
				Particle &p = particles[p_idx];
				if (!p.active) {
					continue;
				}

				float dist = pc.normal.dot(p.position - pc.position);
				if (dist >= 0.0f) {
					continue;
				}

				float penetration = -dist;

				// Unilateral non-penetration constraint: C = dist, grad = normal
				// w = inv_mass, d_lambda = -C / (w + alpha) with alpha = 0 (rigid collision)
				float w = p.inv_mass;
				if (w == 0.0f) continue;
				float d_lambda = penetration / w;
				p.position += d_lambda * w * pc.normal;

				// Friction
				if (pc.friction > 0.0f) {
					Vector3 relative_move = p.position - p.previous_position;
					Vector3 normal_component = relative_move.dot(pc.normal) * pc.normal;
					Vector3 tangent_component = relative_move - normal_component;
					float lateral_dist = tangent_component.length();
					if (lateral_dist > 0.0001f) {
						float friction_reduction = UtilityFunctions::min(pc.friction * penetration, lateral_dist);
						p.position -= tangent_component.normalized() * friction_reduction;
					}
				}
			}
		}

		for (SphereCollider &sc : sphere_colliders) {
			if (!sc.active) {
				continue;
			}

			for (uint32_t p_idx : sc.particles) {
				Particle &p = particles[p_idx];
				if (!p.active) {
					continue;
				}

				Vector3 diff = p.position - sc.position;
				float dist = diff.length();
				float penetration = 0.0f;
				Vector3 normal;

				if (sc.inside) {
					if (dist <= sc.radius) continue;
					normal = dist > 0.0001f ? -diff / dist : Vector3(0, -1, 0);
					penetration = dist - sc.radius;
				} else {
					if (dist >= sc.radius) continue;
					normal = dist > 0.0001f ? diff / dist : Vector3(0, 1, 0);
					penetration = sc.radius - dist;
				}

				// XPBD collision constraint: C = -penetration, grad = normal
				float w = p.inv_mass;
				if (w == 0.0f) continue;
				float d_lambda = penetration / w;
				p.position += d_lambda * w * normal;

				// Friction
				if (sc.friction > 0.0f) {
					Vector3 relative_move = p.position - p.previous_position;
					Vector3 normal_component = relative_move.dot(normal) * normal;
					Vector3 tangent_component = relative_move - normal_component;
					float lateral_dist = tangent_component.length();
					if (lateral_dist > 0.0001f) {
						float friction_reduction = UtilityFunctions::min(sc.friction * penetration, lateral_dist);
						p.position -= tangent_component.normalized() * friction_reduction;
					}
				}
			}
		}

		for (BoxCollider &bc : box_colliders) {
			if (!bc.active) {
				continue;
			}
			for (uint32_t particle_id : bc.particles) {
				Particle &p = particles[particle_id];
				if (!p.active) continue;

				// Transform to local space
				Vector3 local_pos = bc.basis.xform_inv(p.position - bc.position);
				Vector3 half_extents = bc.size * 0.5f;

				// Check if inside
				Vector3 d = local_pos.abs() - half_extents;
				float penetration = 0.0f;
				Vector3 local_normal = Vector3(0, 0, 0);

				if (bc.inside) {
					if (d.x <= 0.0f && d.y <= 0.0f && d.z <= 0.0f) continue;
					// Find the axis with maximum penetration (most outside)
					int axis = 0;
					float max_val = d.x;
					if (d.y > max_val) { axis = 1; max_val = d.y; }
					if (d.z > max_val) { axis = 2; max_val = d.z; }

					penetration = max_val;
					local_normal[axis] = local_pos[axis] > 0 ? -1.0f : 1.0f;
				} else {
					if (d.x >= 0.0f || d.y >= 0.0f || d.z >= 0.0f) continue;
					// Find the axis with minimum penetration (closest face)
					int axis = 0;
					float max_d = d.x;
					if (d.y > max_d) { axis = 1; max_d = d.y; }
					if (d.z > max_d) { axis = 2; max_d = d.z; }

					penetration = -max_d;
					local_normal[axis] = local_pos[axis] > 0 ? 1.0f : -1.0f;
				}

				Vector3 normal = bc.basis.xform(local_normal);

				// XPBD collision constraint
				float w = p.inv_mass;
				if (w == 0.0f) continue;
				float d_lambda = penetration / w;
				p.position += d_lambda * w * normal;

				// Friction
				if (bc.friction > 0.0f) {
					Vector3 relative_move = p.position - p.previous_position;
					Vector3 normal_component = relative_move.dot(normal) * normal;
					Vector3 tangent_component = relative_move - normal_component;
					float lateral_dist = tangent_component.length();
					if (lateral_dist > 0.0001f) {
						float friction_reduction = UtilityFunctions::min(bc.friction * penetration, lateral_dist);
						p.position -= tangent_component.normalized() * friction_reduction;
					}
				}
			}
		}

		// 4. Update velocities and bounciness
		for (uint32_t i = 0; i < particles.size(); i++) {
			Particle &p = particles[i];
			if (!p.active || p.inv_mass == 0.0f) {
				continue;
			}

			p.velocity = (p.position - p.previous_position) * inv_sub_delta;
		}

		// Velocity-based bounciness (post-solver), iterating collider particle lists
		for (PlaneCollider &pc : plane_colliders) {
			if (!pc.active) continue;
			for (uint32_t p_idx : pc.particles) {
				Particle &p = particles[p_idx];
				if (!p.active || p.inv_mass == 0.0f) continue;
				float dist = pc.normal.dot(p.position - pc.position);
				if (dist < 0.01f) {
					float v_dot_n = p.velocity.dot(pc.normal);
					if (v_dot_n < 0.0f) {
						p.velocity -= (1.0f + pc.bounciness) * v_dot_n * pc.normal;
					}
				}
			}
		}

		for (SphereCollider &sc : sphere_colliders) {
			if (!sc.active) continue;
			for (uint32_t p_idx : sc.particles) {
				Particle &p = particles[p_idx];
				if (!p.active || p.inv_mass == 0.0f) continue;
				Vector3 diff = p.position - sc.position;
				float dist = diff.length();
				Vector3 normal;
				bool near_surface = false;

				if (sc.inside) {
					if (dist > sc.radius - 0.01f) {
						near_surface = true;
						normal = dist > 0.0001f ? -diff / dist : Vector3(0, -1, 0);
					}
				} else {
					if (dist < sc.radius + 0.01f) {
						near_surface = true;
						normal = dist > 0.0001f ? diff / dist : Vector3(0, 1, 0);
					}
				}

				if (near_surface) {
					float v_dot_n = p.velocity.dot(normal);
					if (v_dot_n < 0.0f) {
						p.velocity -= (1.0f + sc.bounciness) * v_dot_n * normal;
					}
				}
			}
		}

		for (BoxCollider &bc : box_colliders) {
			if (!bc.active) continue;
			for (uint32_t p_idx : bc.particles) {
				Particle &p = particles[p_idx];
				if (!p.active || p.inv_mass == 0.0f) continue;
				Vector3 local_pos = bc.basis.xform_inv(p.position - bc.position);
				Vector3 half_extents = bc.size * 0.5f;
				Vector3 d = local_pos.abs() - half_extents;
				Vector3 local_normal = Vector3(0, 0, 0);
				bool near_surface = false;

				if (bc.inside) {
					if (d.x > -0.01f || d.y > -0.01f || d.z > -0.01f) {
						near_surface = true;
						int axis = 0;
						float max_val = d.x;
						if (d.y > max_val) { axis = 1; max_val = d.y; }
						if (d.z > max_val) { axis = 2; max_val = d.z; }
						local_normal[axis] = local_pos[axis] > 0 ? -1.0f : 1.0f;
					}
				} else {
					if (d.x < 0.01f && d.y < 0.01f && d.z < 0.01f) {
						near_surface = true;
						int axis = 0;
						float max_d = d.x;
						if (d.y > max_d) { axis = 1; max_d = d.y; }
						if (d.z > max_d) { axis = 2; max_d = d.z; }
						local_normal[axis] = local_pos[axis] > 0 ? 1.0f : -1.0f;
					}
				}

				if (near_surface) {
					Vector3 normal = bc.basis.xform(local_normal);
					float v_dot_n = p.velocity.dot(normal);
					if (v_dot_n < 0.0f) {
						p.velocity -= (1.0f + bc.bounciness) * v_dot_n * normal;
					}
				}
			}
		}
	}
}

void WiggleKitServer::warmup(int p_iterations, const PackedInt32Array &p_particle_ids) {
	if (p_iterations <= 0 || p_particle_ids.is_empty()) {
		return;
	}

	// 1. Back up state
	TightLocalVector<bool> particle_active_backup;
	TightLocalVector<float> particle_inv_mass_backup;
	particle_active_backup.resize(particles.size());
	particle_inv_mass_backup.resize(particles.size());
	for (uint32_t i = 0; i < particles.size(); i++) {
		particle_active_backup[i] = particles[i].active;
		particle_inv_mass_backup[i] = particles[i].inv_mass;
	}

	TightLocalVector<bool> distance_active_backup;
	distance_active_backup.resize(distance_constraints.size());
	for (uint32_t i = 0; i < distance_constraints.size(); i++) {
		distance_active_backup[i] = distance_constraints[i].active;
	}

	TightLocalVector<bool> tetrahedral_active_backup;
	tetrahedral_active_backup.resize(tetrahedral_constraints.size());
	for (uint32_t i = 0; i < tetrahedral_constraints.size(); i++) {
		tetrahedral_active_backup[i] = tetrahedral_constraints[i].active;
	}

	TightLocalVector<bool> plane_active_backup;
	plane_active_backup.resize(plane_colliders.size());
	for (uint32_t i = 0; i < plane_colliders.size(); i++) {
		plane_active_backup[i] = plane_colliders[i].active;
	}

	TightLocalVector<bool> sphere_active_backup;
	sphere_active_backup.resize(sphere_colliders.size());
	for (uint32_t i = 0; i < sphere_colliders.size(); i++) {
		sphere_active_backup[i] = sphere_colliders[i].active;
	}

	TightLocalVector<bool> box_active_backup;
	box_active_backup.resize(box_colliders.size());
	for (uint32_t i = 0; i < box_colliders.size(); i++) {
		box_active_backup[i] = box_colliders[i].active;
	}

	TightLocalVector<bool> plane_sdf_active_backup;
	plane_sdf_active_backup.resize(plane_sdfs.size());
	for (uint32_t i = 0; i < plane_sdfs.size(); i++) {
		plane_sdf_active_backup[i] = plane_sdfs[i].active;
	}

	TightLocalVector<bool> sphere_sdf_active_backup;
	sphere_sdf_active_backup.resize(sphere_sdfs.size());
	for (uint32_t i = 0; i < sphere_sdfs.size(); i++) {
		sphere_sdf_active_backup[i] = sphere_sdfs[i].active;
	}

	TightLocalVector<bool> line_sdf_active_backup;
	line_sdf_active_backup.resize(line_sdfs.size());
	for (uint32_t i = 0; i < line_sdfs.size(); i++) {
		line_sdf_active_backup[i] = line_sdfs[i].active;
	}

	TightLocalVector<bool> box_sdf_active_backup;
	box_sdf_active_backup.resize(box_sdfs.size());
	for (uint32_t i = 0; i < box_sdfs.size(); i++) {
		box_sdf_active_backup[i] = box_sdfs[i].active;
	}

	// 2. Filter participation
	TightLocalVector<bool> is_warmup_particle;
	is_warmup_particle.resize(particles.size());
	for (uint32_t i = 0; i < particles.size(); i++) {
		is_warmup_particle[i] = false;
	}

	for (int i = 0; i < p_particle_ids.size(); i++) {
		uint32_t p_id = (uint32_t)p_particle_ids[i];
		if (p_id < particles.size() && particles[p_id].active) {
			is_warmup_particle[p_id] = true;
		}
	}

	// Deactivate all particles not in warmup
	for (uint32_t i = 0; i < particles.size(); i++) {
		if (!is_warmup_particle[i]) {
			// If it's NOT a warmup particle, it should not MOVE.
			// But it CAN be an obstacle if it was already dynamic.
			if (particles[i].active && particles[i].inv_mass > 0.0f) {
				// Freeze dynamic particles that are not part of warmup
				particles[i].inv_mass = 0.0f;
			}
			// If it was already inactive, it stays inactive.
		}
	}

	// Deactivate constraints/colliders not involving any warmup particle
	for (uint32_t i = 0; i < distance_constraints.size(); i++) {
		if (distance_constraints[i].active) {
			uint32_t pa = distance_constraints[i].particle_a;
			uint32_t pb = distance_constraints[i].particle_b;
			bool involved = false;
			if (pa != INVALID_ID && pa < particles.size() && is_warmup_particle[pa]) involved = true;
			if (pb != INVALID_ID && pb < particles.size() && is_warmup_particle[pb]) involved = true;
			if (!involved) {
				distance_constraints[i].active = false;
			}
		}
	}

	for (uint32_t i = 0; i < tetrahedral_constraints.size(); i++) {
		if (tetrahedral_constraints[i].active) {
			bool involved = false;
			for (int j = 0; j < 4; j++) {
				uint32_t p_id = tetrahedral_constraints[i].particles[j];
				if (p_id != INVALID_ID && is_warmup_particle[p_id]) {
					involved = true;
					break;
				}
			}
			if (!involved) {
				tetrahedral_constraints[i].active = false;
			}
		}
	}

	auto filter_collider = [&](auto &p_collider) {
		if (p_collider.active) {
			bool involved = false;
			for (uint32_t p_id : p_collider.particles) {
				if (is_warmup_particle[p_id]) {
					involved = true;
					break;
				}
			}
			if (!involved) {
				p_collider.active = false;
			}
		}
	};

	for (uint32_t i = 0; i < plane_colliders.size(); i++) filter_collider(plane_colliders[i]);
	for (uint32_t i = 0; i < sphere_colliders.size(); i++) filter_collider(sphere_colliders[i]);
	for (uint32_t i = 0; i < box_colliders.size(); i++) filter_collider(box_colliders[i]);
	for (uint32_t i = 0; i < plane_sdfs.size(); i++) filter_collider(plane_sdfs[i]);
	for (uint32_t i = 0; i < sphere_sdfs.size(); i++) filter_collider(sphere_sdfs[i]);
	for (uint32_t i = 0; i < line_sdfs.size(); i++) filter_collider(line_sdfs[i]);
	for (uint32_t i = 0; i < box_sdfs.size(); i++) filter_collider(box_sdfs[i]);

	// 3. Run simulation
	// Use a fixed delta for stability during warmup
	float warmup_delta = 0.01666f; 
	for (int i = 0; i < p_iterations; i++) {
		iteration(warmup_delta);
	}

	// 4. Restore state
	for (uint32_t i = 0; i < particles.size(); i++) {
		particles[i].active = particle_active_backup[i];
		particles[i].inv_mass = particle_inv_mass_backup[i];
		if (is_warmup_particle[i]) {
			particles[i].velocity = Vector3();
			particles[i].previous_position = particles[i].position;
		}
	}

	for (uint32_t i = 0; i < distance_constraints.size(); i++) {
		distance_constraints[i].active = distance_active_backup[i];
	}

	for (uint32_t i = 0; i < tetrahedral_constraints.size(); i++) {
		tetrahedral_constraints[i].active = tetrahedral_active_backup[i];
	}

	for (uint32_t i = 0; i < plane_colliders.size(); i++) {
		plane_colliders[i].active = plane_active_backup[i];
	}

	for (uint32_t i = 0; i < sphere_colliders.size(); i++) {
		sphere_colliders[i].active = sphere_active_backup[i];
	}

	for (uint32_t i = 0; i < box_colliders.size(); i++) {
		box_colliders[i].active = box_active_backup[i];
	}

	for (uint32_t i = 0; i < plane_sdfs.size(); i++) {
		plane_sdfs[i].active = plane_sdf_active_backup[i];
	}

	for (uint32_t i = 0; i < sphere_sdfs.size(); i++) {
		sphere_sdfs[i].active = sphere_sdf_active_backup[i];
	}

	for (uint32_t i = 0; i < line_sdfs.size(); i++) {
		line_sdfs[i].active = line_sdf_active_backup[i];
	}

	for (uint32_t i = 0; i < box_sdfs.size(); i++) {
		box_sdfs[i].active = box_sdf_active_backup[i];
	}
}
