@tool
class_name WiggleShapeBox extends WiggleShapeBase

@export var Size : Vector3 = Vector3(1, 1, 1):
	set(value):
		Size = value
		update_gizmos()

@export_group("Collider")
@export_range(0.0, 1.0, 0.05) var Bounciness : float = 1
@export_range(0.0, 1.0, 0.05) var Friction : float = 0

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
		_collider_id = WiggleKit.box_collider_create(global_position, global_basis, Size, Bounciness, Friction)
	WiggleKit.box_collider_add_particle(_collider_id, particle_index)

func add_particle_to_sdf_constraint(particle_index:int):
	if Engine.is_editor_hint():
		return
	if _sdf_id == -1:
		_sdf_id = WiggleKit.box_sdf_create(global_position, global_basis, Size, SDF_Distance, SDF_Compliance)
	WiggleKit.box_sdf_add_particle(_sdf_id, particle_index)

func _process(_delta):
	if Engine.is_editor_hint():
		return
	if _collider_id != -1:
		WiggleKit.box_collider_set_position(_collider_id, global_position)
		WiggleKit.box_collider_set_basis(_collider_id, global_basis)
		WiggleKit.box_collider_set_size(_collider_id, Size)
	if _sdf_id != -1:
		WiggleKit.box_sdf_set_position(_sdf_id, global_position)
		WiggleKit.box_sdf_set_basis(_sdf_id, global_basis)
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
