#pragma once

#include "godot_cpp/core/class_db.hpp"
#include "godot_cpp/core/binder_common.hpp"
#include "godot_cpp/templates/local_vector.hpp"
#include "godot_cpp/variant/rid.hpp"
#include "godot_cpp/classes/object.hpp"
#include "godot_cpp/variant/vector3.hpp"
#include "godot_cpp/variant/basis.hpp"

using namespace godot;

class WiggleKitServer : public Object {
	GDCLASS(WiggleKitServer, Object)

	struct Particle {
		Vector3 position;
		Vector3 previous_position;
		Vector3 velocity;
		Vector3 acceleration;
		float inv_mass = 1.0f;
		Vector3 gravity = Vector3(0.0f, -9.81f, 0.0f);
		float damping = 0.0f;
		bool active = false;
	};

	struct DistanceConstraint {
		uint32_t particle_a;
		uint32_t particle_b;
		float rest_length = 0.0f;
		float compliance = 0.0f;
		bool active = false;
	};

	struct PlaneCollider {
		Vector3 position;
		Vector3 normal;
		float bounciness = 0.0f;
		float friction = 0.0f;
		LocalVector<uint32_t> particles;
		bool active = false;
	};

	struct SphereCollider {
		Vector3 position;
		float radius = 1.0f;
		float bounciness = 0.0f;
		float friction = 0.0f;
		LocalVector<uint32_t> particles;
		bool inside = false;
		bool active = false;
	};

	struct BoxCollider {
		Vector3 position;
		Basis basis;
		Vector3 size; // Full width, height, depth
		float bounciness = 0.0f;
		float friction = 0.0f;
		LocalVector<uint32_t> particles;
		bool inside = false;
		bool active = false;
	};

	struct TetrahedralConstraint {
		uint32_t particles[4];
		float rest_volume = 0.0f;
		float compliance = 0.0f;
		bool active = false;
	};

	LocalVector<Particle> particles;
	LocalVector<uint32_t> free_indices;

	LocalVector<DistanceConstraint> distance_constraints;
	LocalVector<uint32_t> free_distance_indices;

	LocalVector<PlaneCollider> plane_colliders;
	LocalVector<uint32_t> free_plane_indices;

	LocalVector<SphereCollider> sphere_colliders;
	LocalVector<uint32_t> free_sphere_indices;

	LocalVector<BoxCollider> box_colliders;
	LocalVector<uint32_t> free_box_indices;

	LocalVector<TetrahedralConstraint> tetrahedral_constraints;
	LocalVector<uint32_t> free_tetrahedral_indices;

	uint32_t substeps = 8;

	static const uint32_t INVALID_ID = 0xFFFFFFFF;
	static WiggleKitServer *singleton;

protected:
	static void _bind_methods();

public:
	static WiggleKitServer *get_singleton();

	WiggleKitServer();
	~WiggleKitServer();

	uint32_t particle_create(const Vector3 &p_position, float p_mass, const Vector3 &p_gravity, float p_damping);
	void particle_free(uint32_t p_id);
	void particle_set_position(uint32_t p_id, const Vector3 &p_position);
	Vector3 particle_get_position(uint32_t p_id) const;
	void particle_set_mass(uint32_t p_id, float p_mass);
	float particle_get_mass(uint32_t p_id) const;
	void particle_set_gravity(uint32_t p_id, const Vector3 &p_gravity);
	Vector3 particle_get_gravity(uint32_t p_id) const;
	void particle_set_damping(uint32_t p_id, float p_damping);
	float particle_get_damping(uint32_t p_id) const;

	uint32_t distance_constraint_create(uint32_t p_particle_a, uint32_t p_particle_b, float p_compliance);
	void distance_constraint_free(uint32_t p_id);

	uint32_t tetrahedral_constraint_create(uint32_t p_particle_a, uint32_t p_particle_b, uint32_t p_particle_c, uint32_t p_particle_d, float p_compliance);
	void tetrahedral_constraint_free(uint32_t p_id);
	void tetrahedral_constraint_set_compliance(uint32_t p_id, float p_compliance);

	uint32_t plane_collider_create(const Vector3 &p_position, const Vector3 &p_normal, float p_bounciness, float p_friction);
	void plane_collider_free(uint32_t p_id);
	void plane_collider_set_position(uint32_t p_id, const Vector3 &p_position);
	void plane_collider_set_normal(uint32_t p_id, const Vector3 &p_normal);
	void plane_collider_add_particle(uint32_t p_id, uint32_t p_particle_id);

	uint32_t sphere_collider_create(const Vector3 &p_position, float p_radius, float p_bounciness, float p_friction, bool p_inside = false);
	void sphere_collider_free(uint32_t p_id);
	void sphere_collider_set_position(uint32_t p_id, const Vector3 &p_position);
	void sphere_collider_set_radius(uint32_t p_id, float p_radius);
	void sphere_collider_add_particle(uint32_t p_id, uint32_t p_particle_id);

	uint32_t box_collider_create(const Vector3 &p_position, const Basis &p_basis, const Vector3 &p_size, float p_bounciness, float p_friction, bool p_inside = false);
	void box_collider_free(uint32_t p_id);
	void box_collider_set_position(uint32_t p_id, const Vector3 &p_position);
	void box_collider_set_basis(uint32_t p_id, const Basis &p_basis);
	void box_collider_set_size(uint32_t p_id, const Vector3 &p_size);
	void box_collider_add_particle(uint32_t p_id, uint32_t p_particle_id);

	void set_substeps(uint32_t p_substeps);
	uint32_t get_substeps() const;

	void iteration(float p_delta);
};
