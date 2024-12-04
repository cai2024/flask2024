

function disableAllButtons() {
    var presetIds = Array.from(document.querySelectorAll('.row > .col-md-12.text-center > .btn') || [])
                         .map(button => button.id)
                         .filter(id => id);
    var buttonIds = [];
    if (typeof available_steps === 'object' && available_steps !== null) {
        Object.keys(available_steps).forEach(function(category) {
            if (Array.isArray(available_steps[category])) {
                available_steps[category].forEach(function(button) {
                    if (button) {
                        buttonIds.push(`btn-${category}-${button}`);
                    }
                });
            }
        });
    }
    var allButtonIds = buttonIds.concat(presetIds);
    allButtonIds.forEach(function(id) {
        var buttonElement = document.getElementById(id);
        if (buttonElement) {
            buttonElement.disabled = true;
        }
    });
}

