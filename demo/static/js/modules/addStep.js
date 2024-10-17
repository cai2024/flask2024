function addStep(type, mode, raw_type) {
    $('#selected-steps').append(`<li class="list-group-item step-item" id="step-${steps.length}">${type}</li>`);
    $('#selected-steps').disableSelection();
    $('#selected-steps').sortable({update: function(event, ui) {
        steps = $('#selected-steps').children().map(function() {
            return steps.find(step => step.id === parseInt(this.id.replace('step-', ''))); }).get();} });
    steps.push({ id: steps.length, type: type, select_option: {}, params:'', mode: mode, raw_type: raw_type, select_index:{}});
}

