@tool
class_name RotationWiggler extends Node3D

# without the factor we would have to set compliance in the range of
# 0.00001 to 0.001...
const COMPLIANCE_FACTOR : float = 0.00001

@export var tip_mass: float = 0.5:
	set(value): 
		if _tip_particle != -1: WiggleKit.particle_set_mass(_tip_particle, value)
		tip_mass = value
@export var tip_distance: float = 0.25
@export var tip_damping : float = 0.5:
	set(value): 
		if _tip_particle != -1: WiggleKit.particle_set_damping(_tip_particle, value)
		tip_damping = value
@export var base_distance: float = 0.75
@export var slackness : float = 0.05:
	set(value): 
		if _tetra_constraint != -1: WiggleKit.tetrahedral_constraint_set_compliance(_tetra_constraint, value * COMPLIANCE_FACTOR)
		slackness = value
@export var gravity : Vector3 = Vector3(0,-10,0):
	set(value): 
		if _tip_particle != -1: WiggleKit.particle_set_gravity(_tip_particle, value)
		gravity = value

var _center_particle: int = -1
var _base_particles: Array[int] = []
var _tip_particle: int = -1

var _tetra_constraint: int = -1
var _dist_constraint: int = -1

# Local positions for static particles
var _local_base_positions: Array[Vector3] = []

func _enter_tree() -> void:
	if Engine.is_editor_hint():
		return
	# the local positions have to be local in the parent space,
	# since our rotation depends on the global positions of the particles.
	var _parent3D := get_parent_node_3d()
	
	# position particle will be kept statically at our global_position
	_center_particle = WiggleKit.particle_create(global_position, 0, Vector3.ZERO, 0)
	
	# 3 base particles forming a tetrahedron base (equilateral triangle in XZ plane)
	var angle_step := TAU / 3.0
	for i in range(3):
		var angle = i * angle_step
		var local_pos := Vector3(cos(angle), 0, sin(angle)) * base_distance
		# we want the starting rotation of ourself to play a role, but we
		# have to work in the local space of our parent, so we do this
		# to_global parent.to_local dance here...
		local_pos = _parent3D.to_local(to_global(local_pos))
		_local_base_positions.append(local_pos)
		var p_id = WiggleKit.particle_create(_parent3D.to_global(local_pos), 0, Vector3.ZERO, 0)
		_base_particles.append(p_id)
		
	# 1 tip particle
	var local_tip_pos = Vector3(0, tip_distance, 0)
	_tip_particle = WiggleKit.particle_create(to_global(local_tip_pos), tip_mass, gravity, tip_damping)
	
	# Constraints
	_tetra_constraint = WiggleKit.tetrahedral_constraint_create(
		_base_particles[0], _base_particles[1], _base_particles[2], _tip_particle, slackness * COMPLIANCE_FACTOR
	)
	_dist_constraint = WiggleKit.distance_constraint_create(_center_particle, _tip_particle, 0.001)

func _exit_tree() -> void:
	if Engine.is_editor_hint():
		return
	if WiggleKit:
		if _tetra_constraint != -1: WiggleKit.tetrahedral_constraint_free(_tetra_constraint)
		if _dist_constraint != -1: WiggleKit.distance_constraint_free(_dist_constraint)
		if _center_particle != -1: WiggleKit.particle_free(_center_particle)
		for p_id in _base_particles:
			WiggleKit.particle_free(p_id)
		if _tip_particle != -1: WiggleKit.particle_free(_tip_particle)

func _physics_process(_delta: float) -> void:
	if Engine.is_editor_hint():
		return
	# first we update our rotation, as that changes the base particle positions
	var tip_pos = WiggleKit.particle_get_position(_tip_particle)
	var base0_pos = WiggleKit.particle_get_position(_base_particles[0])
	
	var up = (tip_pos - global_position).normalized()
	var to_base = (base0_pos - global_position).normalized()
	
	# We want 'up' to be our Y axis.
	# 'to_base' helps define the yaw (rotation around Y).
	
	if up.is_finite() and up.length_squared() > 0.0001:
		var y_axis = up
		# To get a consistent X axis, we cross y_axis with to_base.
		# Since to_base is in the XZ plane locally, this should give us a vector perpendicular to both.
		var x_axis = y_axis.cross(to_base).normalized()
		# Then Z is X cross Y
		var z_axis = x_axis.cross(y_axis).normalized()
		
		if x_axis.length_squared() > 0.0001:
			global_basis = Basis(x_axis, y_axis, z_axis).orthonormalized()
	
	var parent3D := get_parent_node_3d()
	
	# and then we update the static particle positions
	# static center particle
	WiggleKit.particle_set_position(_center_particle, global_position)
	
	# static base particles
	for i in range(3):
		WiggleKit.particle_set_position(_base_particles[i], parent3D.to_global(_local_base_positions[i]))
