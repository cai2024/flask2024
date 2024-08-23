document.getElementById('submitBtn').addEventListener('click', submitSteps);
document.getElementById('unlockBtn').addEventListener('click', function() {
    initialPreset("WGS");
});
document.getElementById('WGSBtn').addEventListener('click', function() {
    initialPreset("WGS");
});

var config_string = "";
var smk_string= "";



async function submitSteps() {
    try {
        const response = await fetch('/submit_form_steps', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json'
            },
            body: JSON.stringify({
                fq_type: fq_type,
                raw_myOptions: raw_myOptions[0],
                steps: steps
            })
        });
        console.log(steps)

        if (response.ok) {
            window.location.href = '/submit_form_code/' + task_id;
        } else {
            const errorData = await response.json();
            console.error('Failed to submit steps:', errorData);
        }
    } catch (error) {
        console.error('Error submitting steps:', error);
    }
}








function disableAllButtons() {
    buttonIds.forEach(function(id) {
        var button = document.getElementById(id);
        if (button) {  // 检查按钮是否存在
            button.disabled = true;  // 禁用按钮
        }
    });
    presetIds.forEach(function(id) {
        var button = document.getElementById(id);
        if (button) {  // 检查按钮是否存在
            button.disabled = true;  // 禁用按钮
        }
    });

}


