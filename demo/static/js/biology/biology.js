$(document).ready(function() {
    initialBiology(form)

});

async function initialBiology(form) {
    try {
        const response = await fetch('/biology/get_preset', {
            method: 'POST', headers: {'Content-Type': 'application/json'},
            body: JSON.stringify({form: form})});
        const data = await response.json();
        data.forEach((step, index) => {
            let html = `<div id="input-${step.step}" class="form-group">
                        <label><strong>${step.step}:</strong></label> `;
            html += `<select id="select-${step.step}" class="form-control mt-2">`;
            step.select.forEach(option => { html += `<option value="${option}">${option}</option>`; });
            html += `</div>`;
            document.getElementById('selected-steps').insertAdjacentHTML('beforeend', html);
        })                        
    } catch (error) { console.error('Error initialPreset:', error);}
}


document.getElementById('confirmBtn').addEventListener('click', confirmSteps);
async function confirmSteps() {
    confirmBtn.disabled = true;  // 禁用确认按钮
    const container = document.getElementById('selected-steps');
    container.querySelectorAll('.form-group').forEach(step => {
        const select = step.querySelector('select');
        const selectedOption = select.options[select.selectedIndex].value;
        steps.push({
            id: steps.length, type: selectedOption,
            mode: form["mode"], raw_type: "raw.fq.gz" ,
        });
    });  
    try {
        const response = await fetch('/batch_code/init_div_steps', {
            method: 'POST', headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ steps :steps })});
        const data = await response.json();
        steps = data
        initOptionsUI(data);
        soft_step = data;
    } catch (error) { console.error('Error:', error); }
    document.querySelectorAll('#selected-steps select').forEach(select => {select.disabled = true;})
    document.getElementById('submitBtn').disabled = false; // 允许提交
}

document.getElementById('submitBtn').addEventListener('click', function() {
    window.location.href = '/query/' + task_id;
    fetch('/submit_smk_job', {
        method: 'POST', headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({ task_id: task_id })
        })
    });



