@tool
extends EditorPlugin

var wiggle_gizmos := WiggleGizmos.new()


func _enable_plugin() -> void:
	add_node_3d_gizmo_plugin(wiggle_gizmos)


func _disable_plugin() -> void:
	remove_node_3d_gizmo_plugin(wiggle_gizmos)


func _enter_tree() -> void:
	# Initialization of the plugin goes here.
	pass


func _exit_tree() -> void:
	# Clean-up of the plugin goes here.
	pass
