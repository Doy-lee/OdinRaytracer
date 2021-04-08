package main

import "core:fmt"
import "core:mem"
import "core:os"

BMPHeader :: struct #packed
{
    // File Header
    file_type: u16, // BM in hex for Windows NT etc
    file_size: u32, // Size of the BMP file in bytes
    reserved: u32, // Application reserved, set to 0 if manually created
    offset_to_pixels: u32, // Byte offset to the pixel data in the image

    // DIB Header: BITMAPINFOHEADER
    dib_header_size: u32, // Must be 40 for BITMAPINFOHEADER
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

    // Pixel Data
    pixels: []byte,
};

BMP_Make :: proc(width, height: i32, bits_per_pixel: u16) -> BMPHeader
{
    result : BMPHeader;
    result.dib_header_size = 40;
    result.width = width;
    result.height = height;
    result.color_planes = 1;
    result.bits_per_pixel = bits_per_pixel;
    result.file_type = 'B' << 8 | 'M' << 0;
    result.pixels = make([]byte, result.width * result.height * cast(i32)result.bits_per_pixel / 8);
    result.file_size = cast(u32)(size_of(result) + len(result.pixels));
    return result;
}

main :: proc()
{
    header := BMP_Make(width = 100, height = 100, bits_per_pixel = 32);
    pixels_rgba := transmute([]u32)header.pixels;

    for _, i in pixels_rgba {
        pixel: []byte = mem.ptr_to_bytes(&pixels_rgba[i], cast(int)header.bits_per_pixel / 8);
        pixel[0] = 255;
        pixel[1] = 0;
        pixel[2] = 0;
        pixel[3] = 255;
    }

    fmt.printf("length {}", len(pixels_rgba));

    OUTPUT_FILE :: "test.bmp";
    bitmap_bytes := mem.ptr_to_bytes(transmute(^u8)&header, cast(int)header.file_size);
    os.write_entire_file(OUTPUT_FILE, bitmap_bytes);
    fmt.printf("File written to {} ({} bytes)", OUTPUT_FILE, bitmap_bytes);
}
