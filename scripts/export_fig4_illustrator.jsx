#target illustrator

(function () {
    var masterFile = new File("D:/桌面/sepsis/文章/返修投稿新主图/Fig4_python_master.pdf");
    var aiFile = new File("D:/桌面/sepsis/文章/返修投稿新主图/Fig4.ai");
    var pdfFile = new File("D:/桌面/sepsis/文章/返修投稿新主图/Fig4.pdf");
    var tifFile = new File("D:/桌面/sepsis/文章/返修投稿新主图/Fig4.tif");
    if (!masterFile.exists) throw new Error("Fig4 master not found");

    app.userInteractionLevel = UserInteractionLevel.DONTDISPLAYALERTS;
    var doc = app.open(masterFile);
    var arial = null;
    try { arial = app.textFonts.getByName("ArialMT"); }
    catch (e1) { try { arial = app.textFonts.getByName("Arial"); } catch (e2) {} }
    if (arial !== null) {
        for (var i = 0; i < doc.textFrames.length; i++) {
            doc.textFrames[i].textRange.characterAttributes.textFont = arial;
        }
    }

    var tifOptions = new ExportOptionsTIFF();
    tifOptions.imageColorSpace = ImageColorSpace.RGB;
    tifOptions.resolution = 300;
    tifOptions.antiAliasing = AntiAliasingMethod.ARTOPTIMIZED;
    tifOptions.lZWCompression = true;
    tifOptions.byteOrder = TIFFByteOrder.IBMPC;
    tifOptions.embedICCProfile = true;
    tifOptions.saveMultipleArtboards = true;
    tifOptions.artboardRange = "1";
    doc.exportFile(tifFile, ExportType.TIFF, tifOptions);
    var artboardTif = new File("D:/桌面/sepsis/文章/返修投稿新主图/Fig4-01.tif");
    if (artboardTif.exists) {
        if (tifFile.exists) tifFile.remove();
        artboardTif.rename("Fig4.tif");
    }

    var pdfOptions = new PDFSaveOptions();
    pdfOptions.compatibility = PDFCompatibility.ACROBAT7;
    pdfOptions.preserveEditability = true;
    pdfOptions.generateThumbnails = true;
    pdfOptions.optimization = true;
    pdfOptions.viewAfterSaving = false;
    doc.saveAs(pdfFile, pdfOptions);

    var aiOptions = new IllustratorSaveOptions();
    aiOptions.pdfCompatible = true;
    aiOptions.embedICCProfile = true;
    aiOptions.compressed = true;
    doc.saveAs(aiFile, aiOptions);
    doc.close(SaveOptions.DONOTSAVECHANGES);

    var transientFiles = [
        masterFile,
        new File("D:/桌面/sepsis/文章/返修投稿新主图/Fig4.png"),
        new File("D:/桌面/sepsis/文章/返修投稿新主图/Fig4.svg"),
        new File("D:/桌面/sepsis/文章/返修投稿新主图/Fig4_preview.png")
    ];
    for (var j = 0; j < transientFiles.length; j++) {
        if (transientFiles[j].exists) transientFiles[j].remove();
    }
    app.userInteractionLevel = UserInteractionLevel.DISPLAYALERTS;
})();
