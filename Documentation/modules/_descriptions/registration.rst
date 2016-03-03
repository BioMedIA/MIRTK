The Registration module of the Medical Image Registration ToolKit (MIRTK)
provides the generic framework used to register images and point sets.
This framework expresses the registration problem as configurable function
minimization problem. The object function, referred to as registration energy
in this context, is put together using the various energy terms. An energy
term can be either an image (dis-)similarity measure, a transformation
constraint, or a point set/surface constraint.
