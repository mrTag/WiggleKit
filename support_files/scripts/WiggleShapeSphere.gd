@tool
class_name WiggleShapeSphere extends WiggleShapeBase

@export var Radius : float = 1:
	set(value):
		Radius = value
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
	# For Node3D, update_gizmos() is not a built-in method.
	# We can use update_configuration_warnings() or just let the editor handle it
	# if the gizmo plugin is watching for changes.
	pass

func add_particle_to_collider(particle_index:int):
	if Engine.is_editor_hint():
		return
	if _collider_id == -1:
		_collider_id = WiggleKit.sphere_collider_create(global_position, Radius, Bounciness, Friction)
	WiggleKit.sphere_collider_add_particle(_collider_id, particle_index)

func add_particle_to_sdf_constraint(particle_index:int):
	if Engine.is_editor_hint():
		return
	if _sdf_id == -1:
		_sdf_id = WiggleKit.sphere_sdf_create(global_position, Radius, SDF_Distance, SDF_Compliance)
	WiggleKit.sphere_sdf_add_particle(_sdf_id, particle_index)

func _process(_delta):
	if Engine.is_editor_hint():
		return
	if _collider_id != -1:
		WiggleKit.sphere_collider_set_position(_collider_id, global_position)
		WiggleKit.sphere_collider_set_radius(_collider_id, Radius)
	if _sdf_id != -1:
		WiggleKit.sphere_sdf_set_position(_sdf_id, global_position)
		WiggleKit.sphere_sdf_set_radius(_sdf_id, Radius)

func _exit_tree() -> void:
	if Engine.is_editor_hint():
		return
	if _collider_id != -1:
		WiggleKit.sphere_collider_free(_collider_id)
		_collider_id = -1
	if _sdf_id != -1:
		WiggleKit.sphere_sdf_free(_sdf_id)
		_sdf_id = -1
