package main

import "core:fmt"
import "core:mem"
import "core:os"

BMPHeader :: struct #packed
{
    file_type: u16, // BM in hex for Windows NT etc
    file_size: u32, // Size of the BMP file in bytes
    reserved01: u16, // Application reserved, set to 0 if manually created
    reserved02: u16, // Application reserved, set to 0 if manually created
    offset_to_pixels: u32, // Byte offset to the pixel data in the image

    // DIB Header
    header_size: u32, // Must be 40 for BITMAPINFOHEADER
    width: i32,
    height: i32,
    color_planes: u16,
    bits_per_pixel: u16,
    compression_method: u32, // Specific magic bytes, or 0 for no compression.
    image_size: u32, // The size of the raw bitmap data. Can be 0 if compression is BI_RGB (0)
    horizontal_resolution: u32, // Horizontal pixels per meter
    vertical_resolution: u32, // Veritcal pixels per meter
    colors_in_color_palette_count: u32, // Or set to 0 to default to 2^n
    important_colors_count: u32, // Or set to 0 when every color is important; generally ignored
};

Dqn_V2i :: [2]int;

BMP_ImageSize :: proc(header: ^BMPHeader) -> Dqn_V2i
{
    result := Dqn_V2i{cast(int)header.width, cast(int)abs(header.height)};
    return result;
}

BMP_Make :: proc(width, height: i32, bits_per_pixel: u16) -> ^BMPHeader
{
    pixels_size_in_bytes := cast(int)(width * height * cast(i32)bits_per_pixel / 8);
    ptr: rawptr = mem.alloc(size = cast(int)size_of(BMPHeader) + pixels_size_in_bytes, alignment = mem.DEFAULT_ALIGNMENT);
    result := transmute(^BMPHeader)ptr;

    result.file_type             = 'M' << 8 | 'B' << 0;
    result.file_size             = cast(u32)(size_of(result^) + pixels_size_in_bytes);
    result.offset_to_pixels      = size_of(result^);
    result.header_size           = 40; // Fixed value for BITMAPINFOHEADER which we support
    result.width                 = width;
    result.height                = -height;
    result.color_planes          = 1;
    result.bits_per_pixel        = bits_per_pixel;
    result.horizontal_resolution = 1;
    result.vertical_resolution   = 1;
    return result;
}

main :: proc()
{
    header := BMP_Make(width = 10, height = 10, bits_per_pixel = 32);
    bytes_per_pixel := cast(int)(header.bits_per_pixel / 8);

    pixels_u8 := mem.ptr_offset(transmute(^u8)header, size_of(header^));
    pixels_size := cast(int)(header.file_size - size_of(header^));

    pixels_slice := mem.slice_ptr(pixels_u8, pixels_size);
    for i := 0; i < len(pixels_slice); i += bytes_per_pixel {
        pixels_slice[i + 0] = 0; // B
        pixels_slice[i + 1] = 0; // G
        pixels_slice[i + 2] = 255; // R
        pixels_slice[i + 3] = 255; // A
        break;
    }

    OUTPUT_FILE :: "test.bmp";
    bitmap_bytes := mem.ptr_to_bytes(transmute(^u8)header, cast(int)header.file_size);
    os.write_entire_file(OUTPUT_FILE, bitmap_bytes);
    fmt.printf("File written to {} (size = {}, bytes = {})", OUTPUT_FILE, header.file_size, bitmap_bytes);
}
