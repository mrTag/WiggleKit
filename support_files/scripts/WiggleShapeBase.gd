@tool
@abstract class_name WiggleShapeBase extends Node3D


# has to be overridden in concrete shape classes
# (add the respective shape constraint to WiggleKit)
@abstract func add_particle_to_collider(particle_index: int)
@abstract func add_particle_to_sdf_constraint(particle_index: int)
