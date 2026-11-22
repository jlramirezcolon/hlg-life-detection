// Print the media-box width and height (points) of a PDF's first page.
import CoreGraphics
import Foundation
let b = CGPDFDocument(URL(fileURLWithPath: CommandLine.arguments[1]) as CFURL)!.page(at: 1)!.getBoxRect(.mediaBox)
print(b.width, b.height)
