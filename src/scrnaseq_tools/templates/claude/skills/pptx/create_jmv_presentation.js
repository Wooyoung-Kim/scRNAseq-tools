const pptxgen = require('pptxgenjs');
const html2pptx = require('./scripts/html2pptx.js');
const path = require('path');

async function createPresentation() {
    const pptx = new pptxgen();
    pptx.layout = 'LAYOUT_16x9';
    pptx.author = 'Claude Code';
    pptx.title = 'JMV_ALI scRNA-seq Analysis Report';
    pptx.subject = 'Multi-Method Cell Type Annotation with Biological Validation';

    const slideDir = '/tmp/claude/-mnt2-kwy-scrna-agent/49993471-81b1-4fb2-b7f2-05757929ae69/scratchpad/pptx_slides';

    console.log('Creating slide 1: Title...');
    await html2pptx(path.join(slideDir, 'slide1.html'), pptx);

    console.log('Creating slide 2: Dataset Overview...');
    await html2pptx(path.join(slideDir, 'slide2.html'), pptx);

    console.log('Creating slide 3: QC Results...');
    await html2pptx(path.join(slideDir, 'slide3.html'), pptx);

    console.log('Creating slide 4: Original Problem...');
    await html2pptx(path.join(slideDir, 'slide4.html'), pptx);

    console.log('Creating slide 5: Multi-Method Approach...');
    await html2pptx(path.join(slideDir, 'slide5.html'), pptx);

    console.log('Creating slide 6: Final Results...');
    await html2pptx(path.join(slideDir, 'slide6.html'), pptx);

    console.log('Creating slide 7: DEG Analysis...');
    await html2pptx(path.join(slideDir, 'slide7.html'), pptx);

    console.log('Creating slide 8: Visualizations...');
    const { slide: slide8, placeholders } = await html2pptx(path.join(slideDir, 'slide8.html'), pptx);

    // Add the multi-method comparison visualization
    if (placeholders.length > 0) {
        const vizPath = '/mnt2/kwy/scrna_agent/results/jmv_ali_epithelial/multi_method_annotation/three_method_comparison_umap.png';
        slide8.addImage({
            path: vizPath,
            ...placeholders[0]
        });
        console.log('  Added visualization image');
    }

    console.log('Creating slide 9: Conclusions...');
    await html2pptx(path.join(slideDir, 'slide9.html'), pptx);

    // Save presentation
    const outputPath = '/mnt2/kwy/scrna_agent/results/jmv_ali_epithelial/JMV_ALI_Analysis_Report.pptx';
    await pptx.writeFile({ fileName: outputPath });
    console.log(`\nPresentation created successfully: ${outputPath}`);
}

createPresentation().catch(err => {
    console.error('Error creating presentation:', err);
    process.exit(1);
});
