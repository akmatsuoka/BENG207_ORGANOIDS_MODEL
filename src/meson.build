cc = meson.get_compiler('c')
m_dep = cc.find_library('m', required : true)

beng297_3_sources = [
  'main.c',
  'beng297_3-application.c',
  'beng297_3-window.c',
]

beng297_3_deps = [
  dependency('gtk4'),
  dependency('libadwaita-1', version: '>= 1.4'),
  m_dep,
]

beng297_3_sources += gnome.compile_resources('beng297_3-resources',
  'beng297_3.gresource.xml',
  c_name: 'beng297_3'
)

executable('beng297_3', beng297_3_sources,
  dependencies: beng297_3_deps,
       install: true,
)
