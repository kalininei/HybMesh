from hybmeshpack import hmscript

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# unit segment
linear_segment = hmscript.create_contour([[0, 0], [1, 0]])

# divide into segments of length 0.03.
# part1 has 33 equal segments
part1 = hmscript.partition_contour(linear_segment, "const", step=0.03)

# partition with refinement towards the center of input line
# segments near end points have length 0.1, at the center - 0.01,
# and between them linear size transition is applied
part2 = hmscript.partition_contour(
    linear_segment, "ref_points",
    step=[0.1, [0, 0], 0.01, [0.5, 0], 0.1, [1.0, 0]])
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

print "partition contour example"
if hmscript.info_contour(part1)['Nedges'] != 33:
    raise Exception
if hmscript.info_contour(part2)['Nedges'] != 26:
    raise Exception
