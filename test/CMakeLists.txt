
add_library(stb_image STATIC stb_image.c stb_image.h stb_image_write.h)

add_executable(compute main.c)
target_link_libraries(compute stb_image distance m)