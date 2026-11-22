// Compose a page: draw each placed PDF (first page, scaled, at x,y in points from the bottom-left), then the overlay.
// Usage: swift place_pdf.swift <overlay.pdf> <out.pdf> [<file.pdf>@x,y,scale ...]
import CoreGraphics
import Foundation

let args = CommandLine.arguments
let overlay = CGPDFDocument(URL(fileURLWithPath: args[1]) as CFURL)!.page(at: 1)!
var box = overlay.getBoxRect(.mediaBox)
let ctx = CGContext(URL(fileURLWithPath: args[2]) as CFURL, mediaBox: &box, nil)!
ctx.beginPDFPage(nil)
for spec in args.dropFirst(3) {
    let parts = spec.components(separatedBy: "@")
    let v = parts[1].split(separator: ",").map { CGFloat(Double($0)!) }
    let page = CGPDFDocument(URL(fileURLWithPath: parts[0]) as CFURL)!.page(at: 1)!
    ctx.saveGState()
    ctx.translateBy(x: v[0], y: v[1])
    ctx.scaleBy(x: v[2], y: v[2])
    ctx.drawPDFPage(page)
    ctx.restoreGState()
}
ctx.drawPDFPage(overlay)
ctx.endPDFPage()
ctx.closePDF()
