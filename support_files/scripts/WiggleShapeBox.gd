@tool
class_name WiggleShapeBox extends WiggleShapeBase

@export var Size : Vector3 = Vector3(1, 1, 1):
	set(value):
		Size = value
		update_gizmos()

@export_group("Collider")
@export_range(0.0, 1.0, 0.05) var Bounciness : float = 1
@export_range(0.0, 1.0, 0.05) var Friction : float = 0
@export_range(0.0, 1.0, 0.05) var VelocityFactor : float = 1

@export_group("SDF Constraint")
@export var SDF_Compliance : float = 0.0
@export var SDF_Distance : float = 0.0

var _collider_id : int = -1
var _sdf_id : int = -1

func update_gizmos():
	pass

func add_particle_to_collider(particle_index:int):
	if Engine.is_editor_hint():
		return
	if _collider_id == -1:
		_collider_id = WiggleKit.box_collider_create(position, basis, Size, Bounciness, Friction)
		WiggleKit.box_collider_set_velocity_factor(_collider_id, VelocityFactor)
	WiggleKit.box_collider_add_particle(_collider_id, particle_index)

func add_particle_to_sdf_constraint(particle_index:int):
	if Engine.is_editor_hint():
		return
	if _sdf_id == -1:
		_sdf_id = WiggleKit.box_sdf_create(position, basis, Size, SDF_Distance, SDF_Compliance)
	WiggleKit.box_sdf_add_particle(_sdf_id, particle_index)

func _physics_process(_delta):
	if Engine.is_editor_hint():
		return
	if _collider_id != -1:
		WiggleKit.box_collider_set_position(_collider_id, position)
		WiggleKit.box_collider_set_basis(_collider_id, basis)
		WiggleKit.box_collider_set_size(_collider_id, Size)
	if _sdf_id != -1:
		WiggleKit.box_sdf_set_position(_sdf_id, position)
		WiggleKit.box_sdf_set_basis(_sdf_id, basis)
		WiggleKit.box_sdf_set_size(_sdf_id, Size)

func _exit_tree() -> void:
	if Engine.is_editor_hint():
		return
	if _collider_id != -1:
		WiggleKit.box_collider_free(_collider_id)
		_collider_id = -1
	if _sdf_id != -1:
		WiggleKit.box_sdf_free(_sdf_id)
		_sdf_id = -1
