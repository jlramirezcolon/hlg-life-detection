// Generic composite merge: draws the original page clipped to the given rectangles (the panels kept unchanged),
// then the overlay page (redrawn panels) on top. Vector throughout.
// Usage: swift merge_keep.swift <original.pdf|-> <overlay.pdf> <out.pdf> [x,y,w,h ...]
//   rectangles are in PDF points with the origin at the bottom-left; "-" as original draws the overlay alone
import CoreGraphics
import Foundation

let args = CommandLine.arguments
let overlay = CGPDFDocument(URL(fileURLWithPath: args[2]) as CFURL)!.page(at: 1)!
var box = overlay.getBoxRect(.mediaBox)
var original: CGPDFPage? = nil
if args[1] != "-" {
    original = CGPDFDocument(URL(fileURLWithPath: args[1]) as CFURL)!.page(at: 1)!
    box = original!.getBoxRect(.mediaBox)
}
let rects = args.dropFirst(4).map { s -> CGRect in
    let v = s.split(separator: ",").map { CGFloat(Double($0)!) }
    return CGRect(x: v[0], y: v[1], width: v[2], height: v[3])
}
let ctx = CGContext(URL(fileURLWithPath: args[3]) as CFURL, mediaBox: &box, nil)!
ctx.beginPDFPage(nil)
if let page = original, !rects.isEmpty {
    ctx.saveGState()
    ctx.clip(to: rects)
    ctx.drawPDFPage(page)
    ctx.restoreGState()
}
let ob = overlay.getBoxRect(.mediaBox)
ctx.saveGState()
ctx.scaleBy(x: box.width / ob.width, y: box.height / ob.height)
ctx.drawPDFPage(overlay)
ctx.restoreGState()
ctx.endPDFPage()
ctx.closePDF()
