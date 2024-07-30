document.getElementById('submitBtn').addEventListener('click', submitSteps);
document.getElementById('unlockBtn').addEventListener('click', function() {
    initialPreset("WGS");
});
document.getElementById('WGSBtn').addEventListener('click', function() {
    initialPreset("WGS");
});




async function submitSteps() {
    console.log('fq_type:', fq_type);
    const stepsData = steps.map(step => ({
        id: step.id,
        type: step.type,
        oninput: step.params,
        select_option: step.select_option,
	select_index: step.select_index,
    }));
    console.log('stepsData:', stepsData);
    try {
        const response = await fetch('/submit_steps', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json'
            },
            body: JSON.stringify({ steps: stepsData ,fq_type:fq_type ,raw_file: raw_myOptions[0] })
        });
        const data = await response.json();         
        if (response.ok) {
            localStorage.setItem('code', data.my_code);
            window.location.href = '/submit_form_code/' + data.task_id; 
          
        } else {
            console.error('Server responded with non-OK status');
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
















