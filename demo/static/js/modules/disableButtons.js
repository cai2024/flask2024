function disableAllButtons() {
    var presetIds = Array.from(document.querySelectorAll('.row > .col-md-12.text-center > .btn')).map(button => button.id);
    var buttonIds = [];
    Object.keys(available_steps).forEach(function(category) {
        available_steps[category].forEach(function(button) {
            buttonIds.push(`btn-${category}-${button}`); });
    });
    var allButtonIds = buttonIds.concat(presetIds);
    allButtonIds.forEach(function(id) {
        if (document.getElementById(id)) {document.getElementById(id).disabled = true;}});
}

