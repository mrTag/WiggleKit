@tool
class_name DanglingWiggler extends Node3D

# without the factor we would have to set compliance in the range of
# 0.00001 to 0.001...
const COMPLIANCE_FACTOR : float = 0.00001

@export var base_mass: float = 0.5:
	set(value): 
		for p_id in _base_particles:
			if p_id != -1: WiggleKit.particle_set_mass(p_id, value/4.0)
		base_mass = value

@export var base_damping : float = 0.5:
	set(value): 
		for p_id in _base_particles:
			if p_id != -1: WiggleKit.particle_set_damping(p_id, value)
		base_damping = value

@export var height: float = 1.0
@export var width: float = 0.5
@export var depth: float = 0.5

@export var gravity : Vector3 = Vector3(0,-10,0):
	set(value): 
		for p_id in _base_particles:
			if p_id != -1: WiggleKit.particle_set_gravity(p_id, value)
		gravity = value

var _position_particle: int = -1
var _base_particles: Array[int] = []

var _tetra_constraint_1: int = -1
var _tetra_constraint_2: int = -1
var _distance_constraints: Array[int] = []

func _enter_tree() -> void:
	if Engine.is_editor_hint():
		return
	
	# position particle will be kept statically at our global_position (tip of pyramid)
	_position_particle = WiggleKit.particle_create(global_position, 999, Vector3.ZERO, 0)
	
	# 4 base particles forming a rectangle in XZ plane, shifted down by 'height'
	# Local coordinates relative to this node
	var base_offsets : Array[Vector3] = [
		Vector3(-width/2, -height, -depth/2),
		Vector3(width/2, -height, -depth/2),
		Vector3(width/2, -height, depth/2),
		Vector3(-width/2, -height, depth/2)
	]
	
	for offset in base_offsets:
		var p_id = WiggleKit.particle_create(to_global(offset), base_mass/4.0, gravity, base_damping)
		_base_particles.append(p_id)
		
	# 2 tetrahedral constraints forming the pyramid
	# Tetra 1: pos, base0, base1, base2
	# Tetra 2: pos, base0, base2, base3
	_tetra_constraint_1 = WiggleKit.tetrahedral_constraint_create(
		_position_particle, _base_particles[0], _base_particles[1], _base_particles[2], 0
	)
	_tetra_constraint_2 = WiggleKit.tetrahedral_constraint_create(
		_position_particle, _base_particles[0], _base_particles[2], _base_particles[3], 0
	)
	
	# Distance constraints between all points to ensure stability
	# 1. Position to all base particles
	for p_id in _base_particles:
		_distance_constraints.append(WiggleKit.distance_constraint_create(_position_particle, p_id, 0))
	
	# 2. Base perimeter
	_distance_constraints.append(WiggleKit.distance_constraint_create(_base_particles[0], _base_particles[1], 0))
	_distance_constraints.append(WiggleKit.distance_constraint_create(_base_particles[1], _base_particles[2], 0))
	_distance_constraints.append(WiggleKit.distance_constraint_create(_base_particles[2], _base_particles[3], 0))
	_distance_constraints.append(WiggleKit.distance_constraint_create(_base_particles[3], _base_particles[0], 0))
	
	# 3. Base diagonals
	_distance_constraints.append(WiggleKit.distance_constraint_create(_base_particles[0], _base_particles[2], 0))
	_distance_constraints.append(WiggleKit.distance_constraint_create(_base_particles[1], _base_particles[3], 0))

func _exit_tree() -> void:
	if Engine.is_editor_hint():
		return
	if WiggleKit:
		if _tetra_constraint_1 != -1: WiggleKit.tetrahedral_constraint_free(_tetra_constraint_1)
		if _tetra_constraint_2 != -1: WiggleKit.tetrahedral_constraint_free(_tetra_constraint_2)
		for c_id in _distance_constraints:
			WiggleKit.distance_constraint_free(c_id)
		_distance_constraints.clear()
		
		if _position_particle != -1: WiggleKit.particle_free(_position_particle)
		for p_id in _base_particles:
			if p_id != -1: WiggleKit.particle_free(p_id)
		_base_particles.clear()

func _physics_process(_delta: float) -> void:
	if Engine.is_editor_hint():
		return
	
	if _base_particles.size() < 4:
		return
		
	var p_pos := global_position
	var b0 := WiggleKit.particle_get_position(_base_particles[0])
	var b1 := WiggleKit.particle_get_position(_base_particles[1])
	var b2 := WiggleKit.particle_get_position(_base_particles[2])
	var b3 := WiggleKit.particle_get_position(_base_particles[3])
	
	var base_center := (b0 + b1 + b2 + b3) / 4.0
	
	# Y axis points from base center to the fixed tip (upwards in pyramid local space)
	var y_axis := (p_pos - base_center).normalized()
	
	# To define yaw, we use the direction between base particles
	# locally base0 to base1 is along X
	var to_side := (b1 - b0).normalized()
	
	if y_axis.is_finite() and y_axis.length_squared() > 0.0001:
		# X axis is perpendicular to Y and our side vector
		var x_axis := y_axis.cross(to_side).cross(y_axis).normalized()
		if x_axis.length_squared() < 0.0001:
			# Fallback if side vector is parallel to Y
			x_axis = Vector3.RIGHT if abs(y_axis.x) < 0.9 else Vector3.FORWARD
			x_axis = x_axis.cross(y_axis).normalized()
			
		var z_axis := x_axis.cross(y_axis).normalized()
		
		global_basis = Basis(x_axis, y_axis, z_axis).orthonormalized()
		
	# Update the static position particle
	WiggleKit.particle_set_position(_position_particle, global_position)
