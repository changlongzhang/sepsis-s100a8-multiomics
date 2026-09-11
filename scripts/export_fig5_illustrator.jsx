#target illustrator
(function(){
 var master=new File("D:/桌面/sepsis/文章/返修投稿新主图/Fig5_python_master.pdf");
 var aiFile=new File("D:/桌面/sepsis/文章/返修投稿新主图/Fig5.ai");
 var pdfFile=new File("D:/桌面/sepsis/文章/返修投稿新主图/Fig5.pdf");
 var tifFile=new File("D:/桌面/sepsis/文章/返修投稿新主图/Fig5.tif");
 app.userInteractionLevel=UserInteractionLevel.DONTDISPLAYALERTS; var doc=app.open(master);
 var arial=null; try{arial=app.textFonts.getByName("ArialMT");}catch(e1){try{arial=app.textFonts.getByName("Arial");}catch(e2){}}
 if(arial!==null){for(var i=0;i<doc.textFrames.length;i++)doc.textFrames[i].textRange.characterAttributes.textFont=arial;}
 var to=new ExportOptionsTIFF(); to.imageColorSpace=ImageColorSpace.RGB; to.resolution=300;
 to.antiAliasing=AntiAliasingMethod.ARTOPTIMIZED; to.lZWCompression=true; to.byteOrder=TIFFByteOrder.IBMPC;
 to.embedICCProfile=true; to.saveMultipleArtboards=true; to.artboardRange="1"; doc.exportFile(tifFile,ExportType.TIFF,to);
 var ab=new File("D:/桌面/sepsis/文章/返修投稿新主图/Fig5-01.tif"); if(ab.exists){if(tifFile.exists)tifFile.remove();ab.rename("Fig5.tif");}
 var po=new PDFSaveOptions(); po.compatibility=PDFCompatibility.ACROBAT7; po.preserveEditability=true; po.optimization=true; po.viewAfterSaving=false; doc.saveAs(pdfFile,po);
 var ao=new IllustratorSaveOptions(); ao.pdfCompatible=true; ao.embedICCProfile=true; ao.compressed=true; doc.saveAs(aiFile,ao); doc.close(SaveOptions.DONOTSAVECHANGES);
 var tmp=[master,new File("D:/桌面/sepsis/文章/返修投稿新主图/Fig5.png"),new File("D:/桌面/sepsis/文章/返修投稿新主图/Fig5.svg"),new File("D:/桌面/sepsis/文章/返修投稿新主图/Fig5_preview.png")];
 for(var j=0;j<tmp.length;j++)if(tmp[j].exists)tmp[j].remove(); app.userInteractionLevel=UserInteractionLevel.DISPLAYALERTS;
})();
