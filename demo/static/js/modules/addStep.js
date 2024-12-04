


function addStep(type, mode, raw_type) {
    let baseOptions = {};
    if (raw_type.endsWith('.fq.gz')) { baseOptions = mode === 'single' ? { fq: [raw_type] }
        : { fq1: [`${raw_type.replace('.fq.gz', '')}_1.fq.gz`], fq2: [`${raw_type.replace('.fq.gz', '')}_2.fq.gz`] };
    } else if (raw_type.endsWith('.bam')) {
        baseOptions = { bam: [raw_type] };
    } else {
        baseOptions = { file: [raw_type] };
    }

    $('#selected-steps').append(`<li class="list-group-item step-item" id="step-${steps.length}">${type}</li>`);
    $('#selected-steps').disableSelection();
    $('#selected-steps').sortable({update: function(event, ui) {
        steps = $('#selected-steps').children().map(function() {
            return steps.find(step => step.id === parseInt(this.id.replace('step-', ''))); }).get();} });
    steps.push({ id: steps.length,type,select_option:{},select_index:{},params:'',mode,baseOptions,selectableOptions:{},rec_option:{},select_index:{}});
}



function addStepOnly(type) {
    $('#selected-steps').append(`<li class="list-group-item step-item" id="step-${steps.length}">${type}</li>`);
    $('#selected-steps').disableSelection();
    $('#selected-steps').sortable({update: function(event, ui) {
        steps = $('#selected-steps').children().map(function() {
            return steps.find(step => step.id === parseInt(this.id.replace('step-', ''))); }).get();} });
}




function initializeStepButtons(buttonSelector, fqTypeSelector, fileTypeSelector) {
    const buttons = document.querySelectorAll(buttonSelector);
    buttons.forEach(button => {
        button.addEventListener('click', function() {
            const buttonType = this.textContent.trim();
            const fqType = document.querySelector(fqTypeSelector).value;
            const fileType = document.querySelector(fileTypeSelector).value;
            addStep(buttonType, fqType, fileType);
        });
    });
}



