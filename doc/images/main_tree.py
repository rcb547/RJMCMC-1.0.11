
import cairo
import random
import math

MM_TO_POINTS = 2.83464567
INCH = int(MM_TO_POINTS * 25.4)
TWOPI = 2.0*math.pi

def draw_labels(ctx, width, top, bottom, thin, med, thick, margin, fontsize, labels, link, step):

    lwidth = width/float(len(labels))
    left = 0.0
    right = lwidth
    cy = (top + bottom)/2.0

    ctx.set_font_size(fontsize)

    i = 0

    for label in labels:

        cw = right - left
        cx = (left + right)/2.0
        
        tw = ctx.text_extents(label)[2]

        ctx.move_to(cx - tw/2.0, cy - step/2.0 + step * float(i % 2))

        ctx.show_text(label)

#        ctx.new_path()
#        ctx.rectangle(left, bottom, right - left, top - bottom)
#        ctx.stroke()

        left = left + lwidth
        right = right + lwidth

        i = i + 1

    tcw = width/float(link[0][0])
    bcw = width/float(link[0][1])
    dh = (bottom - top)/4.0

    ctx.set_line_width(med)

    for l in link[1:]:

        a, b = l

        ctx.new_path()
        ctx.move_to((float(a) + 0.5) * tcw, top - dh)
        ctx.line_to((float(b) + 0.5) * bcw, top + dh)
        if (step > 0.0) and ((b % 2) > 0):
            ctx.line_to((float(b) + 0.5) * bcw, top + 1.5*dh)
        ctx.stroke()

        

    

def draw(ctx, width, height, thin, med, thick, margin, fontsize):

    labels = [['Problem'], 
              ['1D', '2D'],
              ['Single', 'Partitioned', 'Partitioned'],
              ['Regression', 'Forward Model', 'Regression', 'Forward Model', 'Regression', 'Forward Model'],
              ['single1d_regression', 'single1d_direct_regression', 
               'single1d_forwardmodel', 'part1d_regression', 
               'part1d_zero_regression', 'part1d_natural_regression', 
               'part1d_forwardmodel', 'part1d_forwardmodel_natural', 
               'part2d_regression', 'part2d_forwardmodel']]

    links = [[(1,1)],
             [(1,2), (0,0), (0,1)],
             [(2,3), (0,0), (0,1), (1,2)],
             [(3,6), (0,0), (0,1), (1,2), (1,3), (2,4), (2,5)],
             [(6,10), (0,0), (0,1), (1,2), (2,3), (2,4), (2,5), (3,6), (3,7), (4,8), (5,9)]]

    scales = [1.0, 1.0, 0.8, 0.5, 0.5]
    steps = [0.0, 0.0, 0.0, 0.0, 10.0]
    sheight = height/float(len(labels))

    top = 0.0
    bottom = sheight
    scale = 1.0

    for label, scale, link, step in zip(labels, scales, links, steps):
        draw_labels(ctx, width, top, bottom, thin, med, thick, margin, scale * fontsize, label, link, step)
        top = top + sheight
        bottom = bottom + sheight
    
if __name__ == '__main__':

    pdfwidth = 6*INCH
    pdfheight = 3*INCH
    surface = cairo.PDFSurface('main_tree.pdf',
                               pdfwidth,
                               pdfheight)

    ctx = cairo.Context(surface)

    draw(ctx, pdfwidth, pdfheight, 0.1, 0.5, 1.0, 30, 8.0)
    surface.finish()

    pngwidth = 480 * 2
    pngheight = 192 * 2

    surface = cairo.ImageSurface(cairo.FORMAT_RGB24, pngwidth, pngheight)

    ctx = cairo.Context(surface)

    ctx.set_source_rgb(1, 1, 1)
    ctx.new_path()
    ctx.rectangle(0, 0, pngwidth, pngheight)
    ctx.fill()

    ctx.set_source_rgb(0, 0, 0)
    draw(ctx, pngwidth, pngheight, 0.5, 1, 2, 40, 16.0)

    surface.write_to_png('main_tree.png')
