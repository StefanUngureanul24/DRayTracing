TEMPLATE = subdirs

SUBDIRS = QGLViewer OGLRender Raytracing

 # what subproject depends on others
Raytracing.depends = QGLViewer OGLRender
